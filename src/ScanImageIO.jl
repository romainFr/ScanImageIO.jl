module ScanImageIO
using ScanImageTiffReader,FileIO, SharedArrays, Distributed, JSON

## Loads a series of scanImage files, possibly into a binary file (via SharedArrays)
function read_movie_set(files::Array{String,1};binFile=nothing,rois=nothing,channels=nothing,frames=nothing,slices=nothing,volumes=nothing)

    
    @info "Reading metadata" ## First need to check that the movies have compatible sizes
    fullMeta = [ScanImageTiffReader.open(io -> parse_SI_meta(ScanImageTiffReader.metadata(io)),f) for f in files]
    fullSz = [ScanImageTiffReader.open(io -> ScanImageTiffReader.size(io),f) for f in files]
    px = ScanImageTiffReader.open(io -> ScanImageTiffReader.size(io),files[1])
    
    @assert all([f[1:2] == fullSz[1][1:2] for f in fullSz]) "Frames do not have the same size"
    @assert all([frameRate(fm) == frameRate(fullMeta[1]) for fm in fullMeta]) "Not the same frame rate"
    @assert all([hasfastz(fm) == hasfastz(fullMeta[1]) for fm in fullMeta])  "Not all movies are volumetric (or non volumetric)"
    @assert all([volumeRate(fm) == volumeRate(fullMeta[1]) for fm in fullMeta]) "Not the same volume rate"
    @assert all([nrois(fm) == nrois(fullMeta[1]) for fm in fullMeta]) "Not the same ROI structure"
    @assert all([savedChannels(fm) == savedChannels(fullMeta[1]) for fm in fullMeta]) "Not the same channels acquired"
    
    movieSizes = [acquired_channel_frame_slice_volume(fm,fullSz[1]) for fm in fullMeta]
    sz = fullSz[1][1:2]
    nChannels = movieSizes[1][1]
    nFrames = [m[2] for m in movieSizes]
    nAcquiredSlices = [m[3] for m in movieSizes]
    nRealSlices = [nslices(fm) for fm in fullMeta]
    nVolumes = [m[4] for m in movieSizes]
    
    channels = isnothing(channels) ? (1:nChannels) : channels
    frames = [isnothing(frames) ? (1:nFrames[i]) : frames for i in 1:length(movieSizes)]
    slices = [isnothing(slices) ? (1:nRealSlices[i]) : slices for i in 1:length(movieSizes)]
    volumes = [isnothing(volumes) ? (1:nVolumes[i]) : volumes i in 1:length(movieSizes)]
    nRois = nrois(fullMeta[1])
    rois = isnothing(rois) ? (1:nRois) : rois

    nVolumes_Full = sum([length(vol) for vol in volumes])
    
    if nVolumes_Full > length(movieSizes)
        nFrames_Full,nSlices_Full = length(frames[1]),length(slices[1])
        fullDims = 6
    elseif nVolumes_Full == length(movieSizes) ## 1 volume per movie, hence not a volume
        fullDims = 5
        nSlices_Full = sum([length(sl) for sl in slices])
        nFrames_Full = length(frames[1])
        if nSlices_Full == length(movieSizes) ## Repetitions through frames not volumes
            fullDims = 4
            nFrames_Full = sum([length(fr) for fr in frames])
        end
    end
    
    dataLengths = [prod([sz... ,nChannels[i],nFrames[i],nAcquiredSlices[i],nVolumes[i],nRois]) for i in 1:length(fullMeta)]
    
    if hasmRoi(fullMeta[1])
        linesBetweenScanFields = nlinesBetweenFields(fullMeta[1]) 
        roiLines = nlinesPerRoi(fullMeta[1])
        lineOffset = 1
        roisPos = Array{Array{Int64,1},1}(undef,nRois)         
        for i in 1:nRois
            roisPos[i] = lineOffset:(lineOffset+roiLines[i]-1)
            lineOffset = lineOffset + roiLines[i] + linesBetweenScanFields + 1
        end
    else
        roiPos = 1:sz[2]
        roiLines = sz[2]
    end

    dataOut = Array{SharedArray}(undef,length(rois))
    
    i = 1
    if isnothing(binFile)
        for r in rois
            dataOut[i] = SharedArray{px}((roiLines[r],sz[1],length(channels),nFrames_Full,nSlices_Full,nVolumes_Full)[1:fullDims])
            i+=1
        end
    else
        for r in rois
            binFileR = binFile*"_roi_$(r)"
            binFileR = abspath(binFileR)
            dataOut[i] = SharedArray{px}(binFileR,(roiLines[r],sz[1],length(channels),nFrames_Full,nSlices_Full,nVolumes_Full)[1:fullDims])
            i+=1
        end
    end
    
    movingDim = (frames,slices,volumes)[fullD-3]
    offsetMoving = [1;cumsum([length(mD) for mD in movingDim])...]
    @sync begin
        for p in procs(dataOut[1])          
            @async remotecall_wait(write_movie_toshared_native_chunk, p, dataOut,files,nChannels,nFrames,nAcquiredSlices,nVolumes,fullD,channels,frames,slices,volumes,rois,offsetMoving)
        end
    end
    
    dataOut
    
end

## Wrappers for SI metadata
nframes(SImeta) = SImeta["SI"]["hStackManager"]["framesPerSlice"]
nslices(SImeta) = SImeta["SI"]["hStackManager"]["numSlices"]
savedChannels(SImeta) = SImeta["SI"]["hChannels"]["channelSave"]
hasfastz(SImeta) = (SImeta["SI"]["hFastZ"]["hasFastZ"] == 1)
nslicesAcquired(SImeta) =  hasfastz(SImeta) ? SImeta["SI"]["hFastZ"]["numFramesPerVolume"] : nslices(SImeta)
nvolumes(SImeta) = hasfastz(SImeta) ? SImeta["SI"]["hFastZ"]["numVolumes"] : 1
hasmRoi(SImeta) = (SImeta["SI"]["hRoiManager"]["mroiEnable"] == 1)
nrois(SImeta) = hasmRoi(SImeta) ? length(SImeta["RoiGroups"]["imagingRoiGroup"]["rois"]) : 1
volumeRate(SImeta) = SImeta["SI"]["hRoiManager"]["scanVolumeRate"]
frameRate(SImeta) = SImeta["SI"]["hRoiManager"]["scanFrameRate"]
framePeriod(SImeta) = SImeta["SI"]["hRoiManager"]["scanFramePeriod"]
nlinesBetweenFields(SImeta) = round(SImeta["SI"]["hScan2D"]["flytoTimePerScanfield"]/SImeta["SI"]["hRoiManager"]["linePeriod"])
nlinesPerRoi(SImeta) = [SImeta["RoiGroups"]["imagingRoiGroup"]["rois"][i]["scanfields"]["pixelResolutionXY"][2] for i in 1:nrois(SImeta)]

## To deal with cases where the acquisition was aborted
function acquired_channel_frame_slice_volume(SImeta,sz)
    nC = length(savedChannels(SImeta))
    nF = nframes(SImeta)
    nS = nslicesAcquired(SImeta)
    nV = nvolumes(SImeta)
    if sz[3] == nC*nF*nS*nV
        return (nC,nF,nS,nV)
    else
        if sz[3] <= nC*nF
            return (0,0,0,0)
        elseif sz[3] <= nC*nF*nS
            return (nC,nF,div(sz[3],nC*nF),1)
        else
            return (nC,nF,nS,div(sz[3],nC*nF*nS))
        end 
    end
end

## Read a SI movie and return an array of SharedArrays
function write_movie_toshared_chunk(out,files,nChannels,nFrames,nAcquiredSlices,nVolumes,fullD,dataLength,channels,frames,slices,volumes,rois,startPoints,offsets)
   write_movie_toshared_native_chunk(out,files,nChannels,nFrames,nAcquiredSlices,nVolumes,fullD,dataLength,channels,frames,slices,volumes,rois,offsets,myrange(out[1],[0,length(files)]))
end

function write_movie_toshared_chunk(out,files,nChannels,nFrames,nAcquiredSlices,nVolumes,fullD,dataLength,channels,frames,slices,volumes,rois,offsets,filerange)
    for i in filerange
        
        write_movie_to_shared(out,files[i],nChannels,nFrames[i],nAcquiredSlices[i],nVolumes[i],fullD,dataLength[i],channels,frames[i],slices[i],volumes[i],rois,offsets[i])
    end
end

function write_movie_to_shared(out,f::String,sz,nChannels,nFrames,nAcquiredSlices,nVolumes,fullD,dataLength,channels,frames,slices,volumes,rois,offset)
    data = ScanImageTiffReader.open(f) do io
        ScanImageTiffReader.data(io)
    end
    
    for r in rois
        data = reshape(data[1:dataLength],(sz[1],sz[2],nChannels,nFrames,nAcquiredSlices,nVolumes))
        data = permutedims(data,(2,1,3,4,5,6))

        linearOffset = 1 + (offset-1) * stride(out[r],fullD)
        outLength = prod([roisPos[r],sz[1],length(channels),length(frames),length(slices),length(volumes)])
        out[r][linearOffset:(linearOffset+outLength)] = data[roisPos[r],:,channels,frames,slices,volumes]
    end
    out
end

function myrange(q::SharedArray,bigrange)
    idx = indexpids(q)
    if idx == 0 # This worker is not assigned a piece
        return 1:0, 1:0
    end
    nchunks = length(procs(q))
    splits = [round(Int, s) for s in range(bigrange[1],stop=bigrange[2],length=nchunks+1)]
    splits[idx]+1:splits[idx+1]
end




## Reading the metadata header as a Dictionary, in case it's not a JSON header
function makeLineDict(metaLine)
    metaLine = split(metaLine,"= ")
    metaNames = String.(strip.(split(metaLine[1],".")))
    d = Dict{String,Any}(metaNames[end] => parseSIField(metaLine[2]))
    for i in (length(metaNames)-1):-1:1
        d = Dict(metaNames[i] => d)
    end
    d
end

function recursiveMerge(dicts...)
    merge(recursiveMerge,dicts...)
end

function parseSIField(input)
    input = replace(input,"'"=>"")
    out  = try
        eval(Meta.parse(input))
    catch
        input
    end
    out
end

function parse_SI_meta(meta)
    jsonStart = findfirst("{\n",meta)[1]
    if jsonStart == 1
        return JSON.parse(meta)
    else
        rois = JSON.parse(meta[jsonStart:end-1])
        metaP = split(meta[1:(jsonStart-3)],"\n")
        metaP = recursiveMerge([makeLineDict(mL) for mL in metaP]...)
        return merge(metaP,rois)
    end
end

end # module
