module ScanImageIO
using Images,FileIO, ImageMetadata, ImageAxes, ImageMagick, SharedArrays

function scanImage2016Reader(f;shared=false,mapped=false,binFile="imgFile")

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
