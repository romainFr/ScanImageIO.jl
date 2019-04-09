module ScanImageIO
using ScanImageTiffReader,FileIO, SharedArrays, Distributed, JSON

## Loads a series of scanImage files, possibly into a binary file (via SharedArrays)
function scanImage2016Reader(files::Array{String,1};binFile=nothing)

    
    @info "Reading metadata"
    #(pixelsPerLine,linesPerFrame,nSlices,nFrames,resolutionXY,resolutionZ,samplingTime,realSlices) = read_metadata(files[1],read_pos=false)


    
    nFrames = [nFrames]
    
    for f in files[2:end]
        tempMetadata = read_metadata(f,read_pos=false)
        @assert tempMetadata[[1:3;5:7]] == (pixelsPerLine,linesPerFrame,nSlices,resolutionXY,resolutionZ,samplingTime) "Images do not have the same size"
        nFrames = vcat(nFrames,tempMetadata[4])
    end

    ## Finding the longest recording and get the frame positions from this one 
    long_run = findmax(nFrames)[2]
    im_pos = read_metadata(files[long_run],read_pos=true)[9]
 
    nTotalFrames = sum(nFrames)
    if binFile == nothing
        @info "Creating shared Array"
        out = SharedArray{Int16}((linesPerFrame,pixelsPerLine,realSlices,nTotalFrames))
    else
        binFile = abspath(binFile)
        @info "Creating bin file $(binFile)"
        out = SharedArray{Int16}(binFile,(linesPerFrame,pixelsPerLine,realSlices,nTotalFrames))
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

function read_movie(f::String;json=true,rois=nothing,channels=nothing,frames=nothing,slices=nothing,volumes=nothing)
    px,sz,data,metadata = ScanImageTiffReader.open(f) do io
        ScanImageTiffReader.pxtype(io),ScanImageTiffReader.size(io),ScanImageTiffReader.data(io),json ? JSON.parse(ScanImageTiffReader.metadata(io)) : parse_SI_meta(ScanImageTiffReader.metadata(io))
    end 

    fullDims = 6
    channelsAvailable = metadata["SI"]["hChannels"]["channelSave"]
    nChannels = length(channelsAvailable)
    nFrames = metadata["SI"]["hStackManager"]["framesPerSlice"]
    nRealSlices = metadata["SI"]["hStackManager"]["numSlices"]
    if metadata["SI"]["hFastZ"]["hasFastZ"] == 1
        nAcquiredSlices = metadata["SI"]["hFastZ"]["numFramesPerVolume"]
        nVolumes = metadata["SI"]["hFastZ"]["numVolumes"]
    else
        nAcquiredSlices = nRealSlices
        nVolumes = 1
    end

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

    data = reshape(data,(sz[1],sz[2],nChannels,nFrames,nAcquiredSlices,nVolumes))
    data = permutedims(data,(2,1,3,4,5,6))

    if metadata["SI"]["hRoiManager"]["mroiEnable"] == 1
        nRois = length(metadata["RoiGroups"]["imagingRoiGroup"]["rois"])
        linesBetweenScanFields = round(metadata["SI"]["hScan2D"]["flytoTimePerScanfield"]/metadata["SI"]["hRoiManager"]["linePeriod"])
        roiLines = [metadata["RoiGroups"]["imagingRoiGroup"]["rois"][i]["scanfields"]["pixelResolutionXY"][2] for i in 1:nRois]
        lineOffset = 1
        roisPos = Array{Array{Int64,1},1}(undef,nRois)
        rois = isnothing(rois) ? (1:nRois) : rois
        out = Array{Any,1}(undef,length(rois))
        for i in 1:nRois
            roisPos[i] = lineOffset:(lineOffset+roiLines[i]-1)
            lineOffset = lineOffset + roiLines[i] + linesBetweenScanFields + 1
        end
        i = 1
        for r in rois
            out[i] = SharedArray{px}((roiLines[r],sz[1],length(channels),length(frames),length(slices),length(volumes))[1:fullDims])
            out[i][:] = data[roisPos[r],:,channels,frames,slices,volumes]
            i += 1
        end
    else
        out = SharedArray{px}((sz[2],sz[1],length(channels),length(frames),length(slices),length(volumes))[1:fullDims])
        out[:] = data[:,:,channels,frames,slices,volumes]
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
    d = Dict(metaNames[end] => parseSIField(metaLine[2]))
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
