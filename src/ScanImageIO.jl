module ScanImageIO
using Images,FileIO, ImageMetadata, ImageAxes, ImageMagick, SharedArrays, Distributed

function scanImage2016Reader(files::Array{String,1};binFile="imgFile")

    ## TO DO : READ THE HEADER WITHOUT READING THE FULL DATA
    @info "Reading metadata"
    (pixelsPerLine,linesPerFrame,nSlices,nFrames,resolutionXY,resolutionZ,samplingTime,realSlices) = read_metadata(files[1],read_pos=false)
   
    for f in files[2:end]
        tempMetadata = read_metadata(f,read_pos=false)
        @assert tempMetadata[[1:3;5:7]] == (pixelsPerLine,linesPerFrame,nSlices,resolutionXY,resolutionZ,samplingTime) "Images do not have the same size"
        nFrames = vcat(nFrames,tempMetadata[4])
    end

    ## Finding the longest recording and get the frame positions from this one 
    long_run = findmax(nFrames)[2]
    im_pos = read_metadata(files[long_run],read_pos=true)[9]
 
    @info "Creating bin file"
    nTotalFrames = sum(nFrames)
    out = SharedArray{Int16}(binFile,(pixelsPerLine,linesPerFrame,realSlices,nTotalFrames))

    im_length = pixelsPerLine * linesPerFrame * realSlices
    startPoints = vcat(1,cumsum(nFrames*im_length)[1:end-1].-1)

    frame_length = pixelsPerLine * linesPerFrame

    pos_idx = [vcat([i:(i+realSlices-1) for i in 1:nSlices:nF]...) for nF in nFrames]
    im_pos = [im_pos[pI] for pI in pos_idx]
    
    #procsToUse = Distributed.WorkerPool(collect(1:min(max_loads,length(files))))

    #files_split = [round(Int,s) for s in range(0,stop=length(files),length=nchunks+1)]
    #for j in 1:nchunks
    # REDO indexing if more files than proc

    @sync begin
        for p in procs(out)          
            @async remotecall_wait(writeim_toshared_native_chunk, p, out,files,startPoints,im_pos,frame_length)
            #pmap(procsToUse,1:length(files)) do i    
        end
    end
    #end
    ## Returning the SharedArray and a dict of metadata
    (out,Dict("framesPerTrial" => nFrames,"resolutionXY"=>resolutionXY,"resolutionZ"=>resolutionZ,"samplingTime"=>samplingTime))
    
end


function read_metadata(f::String;read_pos=true)
    ## WIP. Just reads the main metadata for now. TOD : MROI situations etc...
    fi = open(f)
    seek(fi,8)
    ifd_header_pos = read(fi,Int)

    @assert read(fi,Int32) == 117637889 "Not a scanimage file" 
    @info "ScanImage TIFF version $(read(fi,Int32))"
    
    nv_frame_data_length = read(fi,Int32)

    roi_group_data_length = read(fi,Int32)

    nv_frame_data = String(read(fi,nv_frame_data_length))

    if read_pos
        ## Collecting where the frames are and their length
        img_pos = []
  
        while ifd_header_pos != 0
            ## Go to the first IFD header
            seek(fi,ifd_header_pos)
            n_ifd = read(fi,Int)
            
            ## Skip to the image location field
            skip(fi,20*6)
            
            @assert read(fi,Int16) == 273 "Failed to find strip offset field"
            skip(fi,10)
            img_pos = vcat(img_pos,read(fi,Int))
        
            ## Finding the next frame
            skip(fi,20*11)
            ifd_header_pos = read(fi,Int)
        end
    end
    
    extracomment = nv_frame_data
    extracomment = split(extracomment,"\n")[1:end-1]
    extracomment = [split(replace(replace(str,"'"=>"\""),"SI."=>""),"= ") for str in extracomment]
    extracomment = Dict(strip(des[1]) => des[2] for des in extracomment)
    
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

    close(fi)
    if read_pos
        return (pixelsPerLine,linesPerFrame,nSlices,nFrames,resolutionXY,resolutionZ,samplingTime,realNSlices,img_pos)
    else
        return (pixelsPerLine,linesPerFrame,nSlices,nFrames,resolutionXY,resolutionZ,samplingTime,realNSlices)
    end
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

function writeim_toshared(out,files,file_split,nFrames,nSlices,startPoints,realSlices)
    filerange = myrange(out,file_split)
    for i in filerange
        img = load(File(format"TIFF",files[i]))
        
        ## Correcting the fact that ImageMagick reads the image as unsigned integers (and return the integer value)
        img = reinterpret(Int16,img)
        img = img[:,:,1:(nFrames[i]*nSlices)]
        
        img = reshape(img,(size(img)[1],size(img)[2],nSlices,nFrames[i]))
        img = img[:,:,1:realSlices,:]
        out[:,:,:,startPoints[i]:(startPoints[i]+nFrames[i]-1)] = img
        @info "Loaded trial $i"
        end
    out
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
   
    fi = open(file)
    for fr in 1:length(im_pos)
        seek(fi,im_pos[fr])
        out[startPoint:(startPoint+imlength-1)] = [read(fi,Int16) for j in 1:imlength]
        startPoint+=imlength
    end
    close(fi)
    @info "Loaded file $(file)"
    out
end


function scanImage2016Reader(f::String;shared=false,mapped=false,binFile="imgFile")

    img = load(File(format"TIFF",f))
    extraprops = magickinfo(f,["comment";"tiff:software"])

    ## Correcting the fact that ImageMagick reads the image as unsigned integers (and return the integer value)
    img = reinterpret(Int16,img)

    comment = extraprops["comment"]
    comment = split(comment,"\n")[1:end-1]
    comment = [split(str,"= ") for str in comment]
    comment = Dict(strip(des[1]) => des[2] for des in comment)
    
    extracomment =extraprops["tiff:software"]
    extracomment = split(extracomment,"\n")[1:end-1]
    extracomment = [split(replace(replace(str,"'","\""),"SI.",""),"= ") for str in extracomment]
    extracomment = Dict(strip(des[1]) => des[2] for des in extracomment)
    
    ## Extract the dimensions/number of frames from the SI comment
    framesAcq = size(img,3)#parse(comment["frameNumberAcquisition"])
    nSlices = 1
    framesPerVolume = parse(extracomment["hFastZ.numFramesPerVolume"])
    if typeof(framesPerVolume) !== Expr
        nSlices = parse(extracomment["hFastZ.numFramesPerVolume"])
    end
    nFrames = div(framesAcq,nSlices)
    if nSlices>1
        samplingTime = (1/parse(extracomment["hRoiManager.scanVolumeRate"]))
    else
        samplingTime = parse(extracomment["hRoiManager.scanFramePeriod"])
    end
   
    img = img[:,:,1:(nFrames*nSlices)]

    img = reshape(img,(size(img)[1],size(img)[2],nSlices,nFrames))
    img = img[:,:,1:(nSlices-parse(extracomment["hFastZ.numDiscardFlybackFrames"])),:]
    #img = Gray.(img)

    resolutionXY = 2*eval(parse(extracomment["hRoiManager.imagingFovUm"]))[2,1]/parse(extracomment["hRoiManager.pixelsPerLine"])
    resolutionZ = parse(extracomment["hStackManager.stackZStepSize"])

    if (shared & !mapped)
        out = SharedArray(img)
    elseif (shared & mapped)
        out = SharedArray{Int16}(binFile,size(img))
        out[:] = img
    end
    
    img = ImageMeta(AxisArray(out,Axis{:x}(range(0,resolutionXY,size(img)[1])),
                              Axis{:y}(range(0,resolutionXY,size(img)[2])),
                              Axis{:z}(range(0,resolutionZ,size(img)[3])),
                              Axis{:time}(range(0,samplingTime,size(img)[4]))),comment=merge(comment,extracomment))
  
    
    img
    
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

