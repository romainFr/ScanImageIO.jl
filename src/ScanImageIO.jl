module ScanImageIO
using Images,FileIO, ImageMetadata, ImageAxes, Unitful, ImageMagick
const um = u"µm"
const s = u"s"

function scanImage2016Reader(f;view=false)

    img = load(File(format"TIFF",f),view=view)
    extraprops = magickinfo(f,["comment";"tiff:software"])

    ## Correcting the fact that ImageMagick reads the image as unsigned integers (and return the integer value)
    #img = reinterpret(Int16,img)

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
    nSlices = parse(extracomment["hFastZ.numFramesPerVolume"])
    nFrames = div(framesAcq,nSlices)
    if nSlices>1
        samplingTime = (1/parse(extracomment["hRoiManager.scanVolumeRate"]))s
    else
        samplingTime = parse(extracomment["hRoiManager.scanFramePeriod"])s
    end
   
    img = img[:,:,1:(nFrames*nSlices)]

    img = reshape(img,(size(img)[1],size(img)[2],nSlices,nFrames))
    img = img[:,:,1:(nSlices-parse(extracomment["hFastZ.numDiscardFlybackFrames"])),:]
    img = Gray.(img)

    resolutionXY = 2um*eval(parse(extracomment["hRoiManager.imagingFovUm"]))[2,1]/parse(extracomment["hRoiManager.pixelsPerLine"])
    resolutionZ = parse(extracomment["hStackManager.stackZStepSize"])um
    
    img = ImageMeta(AxisArray(img,Axis{:x}(range(0um,resolutionXY,size(img)[1])),
                              Axis{:y}(range(0um,resolutionXY,size(img)[2])),
                              Axis{:z}(range(0um,resolutionZ,size(img)[3])),
                              Axis{:time}(range(0s,samplingTime,size(img)[4]))),comment=merge(comment,extracomment))
  
    
    img
    
end

export scanImage2016Reader

end # module
