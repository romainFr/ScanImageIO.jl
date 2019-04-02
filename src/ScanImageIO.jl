module ScanImageIO
using ScanImageTiffReader,FileIO, SharedArrays, Distributed, JSON

## Reading the metadata header as a Dictionary
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
    startPoints = vcat(1,cumsum(nFrames*im_length)[1:end-1].-1)
    
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


function read_metadata(f::String;read_pos=true)
    ## WIP. Just reads the main metadata for now. TOD : MROI situations etc...
    fi = open(f)

    (ifd_header_pos,nv_frame_data_length,roi_group_data_length) = read_SI_header(fi)
   
    extracomment = read_SI_meta(fi,nv_frame_data_length)  

    if read_pos
        ## Collecting where the frames are and their length
        img_pos = []
        while ifd_header_pos != 0
            img_pos_temp,ifd_header_pos = read_ifd_header(fi,ifd_header_pos)
            img_pos = vcat(img_pos,img_pos_temp)
        end
    end

    close(fi)
    
    improps = extract_image_properties(extracomment)
    
    if read_pos
        return (improps...,img_pos)
    else
        return improps
    end
end

function read_ifd_header(fi,ifd_header_pos)
    seek(fi,ifd_header_pos)
    n_ifd = read(fi,Int)
            
    ## Skip to the image location field
    skip(fi,20*6)
    
    @assert read(fi,Int16) == 273 "Failed to find strip offset field"
    skip(fi,10)
    img_pos = read(fi,Int)
    
    ## Finding the next frame
    skip(fi,20*11)
    ifd_header_pos = read(fi,Int)
    return (img_pos,ifd_header_pos)
end

function read_SI_header(fi)
    seek(fi,8)
    ifd_header_pos = read(fi,Int)

    @assert read(fi,Int32) == 117637889 "Not a scanimage file" 
    @info "ScanImage TIFF version $(read(fi,Int32))"
    
    nv_frame_data_length = read(fi,Int32)

    roi_group_data_length = read(fi,Int32)
    return (ifd_header_pos,nv_frame_data_length,roi_group_data_length) 
end

function read_SI_meta(fi,header_length)
    parse_SI_meta(String(read(fi,header_length)))
end



function extract_image_properties(extracomment)
    framesPerVolume = Meta.parse(extracomment["hFastZ.numFramesPerVolume"])
    if typeof(framesPerVolume) !== Expr
        nSlices = Meta.parse(extracomment["hFastZ.numFramesPerVolume"])
        nFrames = Meta.parse(extracomment["hFastZ.numVolumes"])
        samplingTime = (1/Meta.parse(extracomment["hRoiManager.scanVolumeRate"]))
    else
        nFrames = Meta.parse(extracomment["hStackManager.framesPerSlice"])
        nSlices = 1
        samplingTime = Meta.parse(extracomment["hRoiManager.scanFramePeriod"])
    end
    pixelsPerLine = Meta.parse(extracomment["hRoiManager.pixelsPerLine"])
    linesPerFrame = Meta.parse(extracomment["hRoiManager.linesPerFrame"])
    realNSlices = nSlices - Meta.parse(extracomment["hFastZ.numDiscardFlybackFrames"])
    
    resolutionXY = 2*eval(Meta.parse(extracomment["hRoiManager.imagingFovUm"]))[2,1]/Meta.parse(extracomment["hRoiManager.pixelsPerLine"])
    resolutionZ = Meta.parse(extracomment["hStackManager.stackZStepSize"])

    return (pixelsPerLine,linesPerFrame,nSlices,nFrames,resolutionXY,resolutionZ,samplingTime,realNSlices)
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

function scanImage5Reader(f;objMag=40,shared=false,mapped=false,binFile="imgFile")

    img = load(File(format"TIFF",f))
    extraprops = magickinfo(f,["comment"])

    ## Correcting the fact that ImageMagick reads the image as unsigned integers (and return the integer value)
    img = reinterpret(Int16,img)

    comment = extraprops["comment"]
    comment = split(comment,"\n")[1:end-1]
    comment = [split(replace(replace(str,"'","\""),"scanimage.SI5.","")," = ") for str in comment]
    comment = Dict(strip(des[1]) => des[2] for des in comment)
 
    ## Extract the dimensions/number of frames from the SI comment
    framesAcq = size(img,3)#parse(comment["frameNumberAcquisition"])
    nSlices = nSlices =  parse(comment["stackNumSlices"])
    nFrames = div(framesAcq,nSlices)
    if nSlices>1
        samplingTime = samplingTime = parse(comment["fastZPeriod"])
    else
        samplingTime = samplingTime = parse(comment["scanFramePeriod"])
    end
   
    img = img[:,:,1:(nFrames*nSlices)]

    img = reshape(img,(size(img)[1],size(img)[2],nSlices,nFrames))
    img = img[:,:,1:(nSlices-parse(comment["fastZNumDiscardFrames"])),:]
    #img = Gray.(img)
    
    resolutionXY = ((tan(10*2*pi/360)*45000/objMag)/(parse(comment["pixelsPerLine"])*parse(comment["zoomFactor"])))
    resolutionZ = (parse(comment["stackZStepSize"])/10)
   
    if shared
        img = SharedArray(img)
    elseif mapped
        img = SharedArray(binFile,img)
    end
    
    img = ImageMeta(AxisArray(img,Axis{:x}(range(0,resolutionXY,size(img)[1])),
                              Axis{:y}(range(0,resolutionXY,size(img)[2])),
                              Axis{:z}(range(0,resolutionZ,size(img)[3])),
                              Axis{:time}(range(0,samplingTime,size(img)[4]))),comment=comment)
  
    
    img
    
end
 

export scanImage2016Reader, scanImage5Reader

end # module

