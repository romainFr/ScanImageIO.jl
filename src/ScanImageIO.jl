module ScanImageIO
using Images,FileIO


function scanImage2016Reader(f)

    img = load(File(format"TIFF",f),extraprop =["comment";"tiff:software"])

    ## Correcting the fact that ImageMagick reads the image as unsigned integers (and return the integer value)
    img = reinterpret(Int16,img)

    comment = img["comment"]
    comment = split(comment,"\n")[1:end-1]
    comment = [split(str,"= ") for str in comment]
    comment = [strip(des[1]) => des[2] for des in comment]
    
    extracomment =img["tiff:software"]
    extracomment = split(extracomment,"\n")[1:end-1]
    extracomment = [split(replace(replace(str,"'","\""),"SI.",""),"= ") for str in extracomment]
    extracomment = [strip(des[1]) => des[2] for des in extracomment]
    
    ## Extract the dimensions/number of frames from the SI comment
    framesAcq = parse(comment["frameNumberAcquisition"])
    nSlices = parse(extracomment["hFastZ.numFramesPerVolume"])
    nFrames = div(framesAcq,nSlices)
    if nSlices>1
        samplingTime = 1/parse(extracomment["hRoiManager.scanVolumeRate"])     
    else
        samplingTime = parse(extracomment["hRoiManager.scanFramePeriod"])
    end
   
    img = img[:,:,1:(nFrames*nSlices)]

    img = reshape(img,(size(img)[1],size(img)[2],nSlices,nFrames))
    img = img[:,:,1:(nSlices-parse(extracomment["hFastZ.numDiscardFlybackFrames"])),:]
    img = grayim(Image(img))
    img["comment"] = merge(comment,extracomment)
    ## Input the correct properties
    img["timedim"] = 4
    img["spatialorder"]=["x","y","z"]
    #resolutionXY = (tan(10*2*pi/360)*45000/objMag)/(parse(extracomment["hRoiManager.pixelsPerLine"])*parse(extracomment["hRoiManager.scanZoomFactor"]))
    resolutionXY = 2*eval(parse(extracomment["hRoiManager.imagingFovUm"]))[2,1]/parse(extracomment["hRoiManager.pixelsPerLine"])
    resolutionZ = parse(extracomment["hStackManager.stackZStepSize"])
    
    img["pixelspacing"] = [resolutionXY,resolutionXY,resolutionZ]
    img["period"]=samplingTime
    img["colorspace"]="Gray"
    img
    
end

export scanImage2016Reader

end # module
