module ScanImageIO
using ScanImageTiffReader,FileIO, SharedArrays, Distributed, JSON

## Loads a series of scanImage files, possibly into a binary file (via SharedArrays)
function read_movie_set(files::Array{String,1};binFile=nothing,json=true)

    
    @info "Reading metadata" ## First need to check that the movies have compatible sizes
    fullMeta = [ScanImageTiffReader.open(io -> json ? JSON.parse(ScanImageTiffReader.metadata(io)) : parse_SI_meta(ScanImageTiffReader.metadata(io)),f)]
    fullSz = [ScanImageTiffReader.open(io -> ScanImageTiffReader.size(io),f)]
    
    @assert all([f[1:2] == fullSz[1][1:2] for f in fullSz]) "Frames do not have the same size"
    @assert all([frameRate(fm) == frameRate(fullMeta[1]) for fm in fullMeta]) "Not the same frame rate"
    @assert all([hasfastz(fm) == hasfastz(fullMeta[1]) for fm in fullMeta])  "Not all movies are volumetric (or non volumetric)"
    @assert all([volumeRate(fm) == volumeRate(fullMeta[1]) for fm in fullMeta]) "Not the same volume rate"
    @assert all([nrois(fm) == nrois(fullMeta[1]) for fm in fullMeta]) "Not the same ROI structure"
    
    movieSizes = [acquired_channel_frame_slice_volume(fm,fullSz[1]) for fm in fullMeta]

    if hasfastz(fullMeta[1])
        nTotalVolumes = sum([mS[4] for mS in movieSizes])
        if binFile == nothing
            @info "Creating shared Array"
            out = SharedArray{Int16}((fullSz[1][2],fullSz[1][1],nslices(fullMeta[1]),nTotalVolumes))
        else
            binFile = abspath(binFile)
            @info "Creating bin file $(binFile)"
            out = SharedArray{Int16}(binFile,(fullSz[1][2],fullSz[1][1],nslices(fullMeta[1]),nTotalVolumes))
        end
    end
        
    im_length = pixelsPerLine * linesPerFrame * realSlices
    
    frame_length = pixelsPerLine * linesPerFrame

    pos_idx = [vcat([i:(i+realSlices-1) for i in 1:nSlices:nF]...) for nF in nFrames*nSlices]
    im_pos = [im_pos[pI] for pI in pos_idx]
    
    @sync begin
        for p in procs(out)          
            @async remotecall_wait(writeim_toshared_native_chunk, p, out,files,startPoints,im_pos,frame_length)
        end
    end
 
    ## Returning the SharedArray and a dict of metadata
    (out,Dict("framesPerTrial" => nFrames,"resolutionXY"=>resolutionXY,"resolutionZ"=>resolutionZ,"samplingTime"=>samplingTime))
    
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
function read_movie(f::String;json=true,rois=nothing,channels=nothing,frames=nothing,slices=nothing,volumes=nothing,binFile=nothing)
    px,sz,data,metadata = ScanImageTiffReader.open(f) do io
        ScanImageTiffReader.pxtype(io),ScanImageTiffReader.size(io),ScanImageTiffReader.data(io),json ? JSON.parse(ScanImageTiffReader.metadata(io)) : parse_SI_meta(ScanImageTiffReader.metadata(io))
    end 

    fullDims = 6
    channelsAvailable = savedChannels(metadata)
    nChannels,nFrames,nAcquiredSlices,nVolumes = acquired_channel_frame_slice_volume(metadata,sz)
    dataLength = sz[1]*sz[2]*nChannels*nFrames*nAcquiredSlices*nVolumes   
    nRealSlices = nslices(metadata)


    if nVolumes == 1
        fullDims =5
        if nRealSlices == 1
            fullDims = 4
        end
    end

    channels = isnothing(channels) ? (1:nChannels) : channels
    frames = isnothing(frames) ? (1:nFrames) : frames
    slices = isnothing(slices) ? (1:nRealSlices) : slices
    volumes = isnothing(volumes) ? (1:nVolumes) : volumes

    data = reshape(data[1:dataLength],(sz[1],sz[2],nChannels,nFrames,nAcquiredSlices,nVolumes))
    data = permutedims(data,(2,1,3,4,5,6))

    if hasmRoi(metadata)
        nRois = nrois(metadata)
        linesBetweenScanFields = nlinesBetweenFields(metadata) 
        roiLines = nlinesPerRoi(metadata)
        lineOffset = 1
        roisPos = Array{Array{Int64,1},1}(undef,nRois)
        rois = isnothing(rois) ? (1:nRois) : rois
        dataOut = Array{Any,1}(undef,length(rois))
        for i in 1:nRois
            roisPos[i] = lineOffset:(lineOffset+roiLines[i]-1)
            lineOffset = lineOffset + roiLines[i] + linesBetweenScanFields + 1
        end
        i = 1
        if isnothing(binFile)
            for r in rois
                dataOut[i] = SharedArray{px}((roiLines[r],sz[1],length(channels),length(frames),length(slices),length(volumes))[1:fullDims])
                dataOut[i][:] = data[roisPos[r],:,channels,frames,slices,volumes]
                i += 1
            end
        else
            for r in rois
                binFileR = binFile*"_roi_$(r)"
                binFileR = abspath(binFileR)
                dataOut[i] = SharedArray{px}(binFileR,(roiLines[r],sz[1],length(channels),length(frames),length(slices),length(volumes))[1:fullDims])
                dataOut[i][:] = data[roisPos[r],:,channels,frames,slices,volumes]
                i += 1
            end
        end
    else
        dataOut = Array{Any,1}(undef,1)
        if isnothing(binFile)
            dataOut[1] = SharedArray{px}((sz[2],sz[1],length(channels),length(frames),length(slices),length(volumes))[1:fullDims])
            dataOut[1][:] = data[:,:,channels,frames,slices,volumes]
        else
            binFile = abspath(binFile)
            dataOut[1] = SharedArray{px}(binFile,(sz[2],sz[1],length(channels),length(frames),length(slices),length(volumes))[1:fullDims])
        end
    end
    dataOut
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

function writeim_toshared_native_chunk(out,files,startPoints,im_pos,imlength)
   writeim_toshared_native_chunk(out,files,startPoints,im_pos,imlength,myrange(out,[0,length(files)]))
end

function writeim_toshared_native_chunk(out,files,startPoints,im_pos,imlength,filerange)
    for i in filerange
        writeim_toshared_native(out,files[i],startPoints[i],im_pos[i],imlength)
    end
end

function writeim_toshared_native(out,file,startPoint,im_pos,imlength)

    imsz = (size(out,2),size(out,1))
    
    try
        fi = open(file)
        sP = startPoint
        for fr in 1:length(im_pos)
            seek(fi,im_pos[fr])
            out[sP:(sP+imlength-1)] = permutedims(reshape([read(fi,Int16) for j in 1:imlength],imsz),(2,1))
            sP+=imlength
        end
        @info "Loaded file $(file)"
        close(fi)
        return out
    catch
        @info "Need to read the image positions for file $(file)"
        metad = read_metadata(file,read_pos=true)
        realSlices = metad[8]
        im_pos = metad[9]
        nSlices = metad[3]
        nFrames =metad[4]
        pos_idx = vcat([i:(i+realSlices-1) for i in 1:nSlices:nFrames]...)
        im_pos = im_pos[pos_idx]
        fi = open(file)
        for fr in 1:length(im_pos)
            seek(fi,im_pos[fr])
            out[startPoint:(startPoint+imlength-1)] = permutedims(reshape([read(fi,Int16) for j in 1:imlength],imsz),(2,1))
            startPoint+=imlength
        end
        @info "Loaded file $(file)"
        close(fi)
        return out
    end
   
   
end

## Reads a single file.
function scanImage2016Reader(file::String;binFile=nothing)
    scanImage2016Reader([file];binFile=binFile)
end

export scanImage2016Reader, scanImage5Reader


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
    rois = JSON.parse(meta[jsonStart:end-1])
    metaP = split(meta[1:(jsonStart-3)],"\n")
    metaP = recursiveMerge([makeLineDict(mL) for mL in metaP]...)
    merge(metaP,rois)
end

end # module
