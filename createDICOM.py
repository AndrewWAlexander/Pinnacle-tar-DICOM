###################################################################################################
#                                    Program:  createdcm.py                                       #
#    This program was written to read in pinnacle text files and convert them into DICOM files    #
#    The files required: Patient (text file under patient folder),                                #
#                        plan.Points (file under Patient folder->plan label folder),              #
#                        plan.roi (file under Patient folder -> plan label folder),               #
#                        plan.Trial (file under Patient folder -> plan label folder),             #
#                        plan.Trial.Binary.xxx (xxx represents number of file, several needed)    #
#                        ImageSet_%s.ImageInfo (text file under patient folder)                   #
#                        ImageSet_%s.header                                                       #
#                        Either folder with images or ImageSet_%s file
#						 Code Authors are: Colleen Henschel,  Andrew Alexander                    #
###################################################################################################


####################################################################################################################################################
#   import libraries below
####################################################################################################################################################
from __future__ import print_function

import time #used for getting current date and time for file
import re #used for isolated values from strings
import sys
import os.path
import pydicom
import numpy as np
from pydicom.dataset import Dataset, FileDataset
from pydicom.sequence import Sequence
from pydicom.filebase import DicomFile
import pydicom.uid
import os
import struct
from random import randint
from PIL import Image

####################################################################################################################################################
#  Global Variables
####################################################################################################################################################
ROI_COUNT = 0 #This value will represent the ROI that I'm currently looking at in file, will be incremented for each roi
SeriesUID = "NA"
StudyInstanceUID = "NA"
FrameUID = "NA"
ClassUID = "NA"
patientname = ""
dob = ""
pid = ""
imageslice = []
imageuid = []
Colors = [['255','0','0'],['255','20','147'],['0','0','255'],['0','255','0'],['125','38','205'],['255','255','0'],['255','140','0'],['0','100','0'],['0','191','255'],['255','192','203'],['72','209','204'],['139','69','19'],['255','193','37'],['221','160','221'],['107','142','35'],['142','35','35'],['245','204','176'],['191','239','255'],['139','28','98'],['255','99','71']] #red, pink, blue, green, purple, yellow, orange, dark green, sky blue, light pink, Turquois, brown, gold,lightpurple, olive, brick, peach?, light blue, maroon, tomato  
patient_sex = ""
study_date = ""
study_time = ""
model = ""
physician = ""
sid = ""
isocenter = []
ctcenter = []
descrip = ""
plancount = 0
plannamelist = []
planids = []
randval = randint(0,999)
currentdate = time.strftime("%Y%m%d")
currenttime = time.strftime("%H%M%S")
doserefpt = []
patient_position = ""
xshift = 0
yshift = 0
zshift = 0
lname = ""
fname = ""
patientfolder = ""
structsopinstuid = ''
structseriesinstuid = ''
plansopinstuid = ''
planseriesinstuid = ''
doseseriesuid = ''
doseinstuid = ''
planfilename = ''
dosexdim = 0
doseydim = 0
dosezdim = 0
doseoriginx = ""
doseoriginy = ""
doseoriginz = ""
beamdosefiles = []
pixspacingx = ""
pixspacingy = ""
pixspacingz = ""
posrefind = ""
image_orientation = []
imagesetnumber = ""
point_names = []
point_values = []
numfracs = ""
flag_nobinaryfile = False
flag_noimages = False
GTransferSyntaxUID='1.2.840.10008.1.2'
no_setup_file = False
no_beams = False
gImplementationClassUID='1.2.826.0.1.3680043.8.498.75006884747854523615841001'
Manufacturer="Pinnacle Philips"
PDD6MV = 0.6683
PDD10MV = 0.6683 # Also temporary, need to get actual PDD value
PDD15MV = 0.7658
PDD16MV = 0.7658 ## THis is temporary, this value is not correct just using as place holder for now
softwarev = ""
slicethick = 0
x_dim = 0
y_dim = 0
z_dim = 0
xpixdim = 0
ypixdim = 0
#listofversions = []

####################################################################################################################################################
# Function: main
# This main function is what should be called to run program
# The name of the patient folder should be passed into the function for it to be run. This folder should be stored under Inputf directory (line 102)
####################################################################################################################################################
def main(temppatientfolder,inputfolder,outputfolder):
    global ROI_COUNT
    global structfilename
    global SeriesUID
    global StudyInstanceUID
    global FrameUID
    global ClassUID
    global patientname
    global dob
    global pid
    global imageslice
    global imageuid
    global Colors  
    global patient_sex
    global study_date
    global study_time
    global model
    global physician
    global sid
    global isocenter
    global ctcenter
    global descrip
    global plancount
    global plannamelist
    global planids
    global randval
    global currentdate
    global currenttime
    global doserefpt
    global patient_position
    global xshift
    global yshift
    global lname
    global fname
    global Inputf
    global Outputf
    global patientfolder
    global structsopinstuid
    global structseriesinstuid
    global plansopinstuid
    global planseriesinstuid
    global doseinstuid
    global doseseriesuid
    global planfilename
    global flag_noimages
    global no_setup_file
    global no_beams
    global softwarev
    
    initglobalvars()  # First step of the main function is to call the initglobalvars variable to reset everything in case this function is being used in a loop. (see allpatientloop.py) 
    #print("Input Patient Folder:")
    #patientfolder = raw_input("> ")
    
    patientfolder = temppatientfolder 
    Inputf = inputfolder
    Outputf = outputfolder


    if not os.path.exists(Outputf+"%s"%(patientfolder)):
        os.makedirs(Outputf+"%s"%(patientfolder)) #Create folder for exported DICOM files if it does not already exist

    sys.stdout = open(Outputf + "%s/log.txt"%patientfolder, "w")


    print("Pinnacle tar folder path: " + Inputf)
    print("Current Patient: " + patientfolder)
    print("Output location: " + Outputf) 
    
    structsopinstuid = pydicom.uid.generate_uid() 
    structds = createstructds() #creating dataset for structure file
    for j in range(0, 5000):
        morewastingtime = j
    structseriesinstuid = pydicom.uid.generate_uid()
    structds.ReferencedStudySequence = Sequence()
    
    structds = initds(structds)

    structds = readpatientinfo(structds)
    readImageInfo() #Gets UID information for image files 
    structds = initds(structds) #initializes values like uids, creation time, manufacturer, values that are not patient dependent
    
    for i in range(0, 5000):
        timewaster = i
    plansopinstuid = pydicom.uid.generate_uid()
    convertimages() #This function makes the image files usable (matches patient info that will go into other DICOM files). If image files do not exist this function calls createimagefiles function
    if flag_noimages:
        return
    planseriesinstuid = pydicom.uid.generate_uid()

    patient_position = getpatientsetup("Plan_%s"%planids[0])
    if no_setup_file == True:
        return

    structds.ReferencedFrameOfReferenceSequence = Sequence()
    ReferencedFrameofReference1 = Dataset()
    structds.ReferencedFrameOfReferenceSequence.append(ReferencedFrameofReference1)
    structds.ReferencedFrameOfReferenceSequence[0].FrameofReferenceUID = FrameUID
    structds.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence = Sequence()
    RTReferencedStudy1 = Dataset()
    structds.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence.append(RTReferencedStudy1)
    structds.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].ReferencedSOPClassUID = '1.2.840.10008.3.1.2.3.2'
    structds.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].ReferencedSOPInstanceUID = StudyInstanceUID
    structds.StudyInstanceUID = StudyInstanceUID
    structds.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].RTReferencedSeriesSequence = Sequence()
    RTReferencedSeries1 = Dataset()
    structds.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].RTReferencedSeriesSequence.append(RTReferencedSeries1)
    structds.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].RTReferencedSeriesSequence[0].SeriesInstanceUID = SeriesUID
    structds.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].RTReferencedSeriesSequence[0].ContourImageSequence = Sequence()
    for i, value in enumerate(imageuid,1): 
        exec("ContourImage%d = Dataset()"% i)
        exec("structds.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].RTReferencedSeriesSequence[0].ContourImageSequence.append(ContourImage%d)"%i)
        structds.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].RTReferencedSeriesSequence[0].ContourImageSequence[i-1].ReferencedSOPClassUID = '1.2.840.10008.5.1.4.1.1.1'
        structds.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].RTReferencedSeriesSequence[0].ContourImageSequence[i-1].ReferencedSOPInstanceUID = imageuid[i-1]

    doseinstuid = pydicom.uid.generate_uid()

    structds.ROIContourSequence = Sequence()
    structds.StructureSetROISequence = Sequence()
    structds.RTROIObservationsSequence = Sequence()
    
    if softwarev != "Pinnacle 9.0": # If pinnacle software is version 9.0 the shifts are not needed so this function can be skipped and the values for the shifts will still be set to zero
        getstructshift()
    structds = readpoints(structds, "Plan_%s"%planids[0])
    structds = readroi(structds, "Plan_%s"%planids[0])

    structds.ApprovalStatus = 'UNAPPROVED' #find out where to get if its been approved or not
    # Set the transfer syntax
    structds.is_little_endian = True
    structds.is_implicit_VR = True
    #structfilepath=outputfolder + patientfolder + "/" + structfilename
    #structds.save_as("structfilepath")
    print("Structure file being saved\n")
    structds.save_as(Outputf + "/%s/%s"%(patientfolder, structfilename))
    doseseriesuid = pydicom.uid.generate_uid()
    print("creating plan data structures \n")
    
    #############################################################################################
    # loop below creates plan files for each plan in directory (based on what is in the Patient file)
    for i in range(0, plancount): 
        
        planame = plannamelist[i]
        plandirect = "Plan_" + planids[i]
        exec("plands_%s = createplands(i)"%planids[i])
        exec("plands_%s = planinit(plands_%s, planame, plandirect, i)"%(planids[i], planids[i]))
        exec("plands_%s = readtrial(plands_%s, plandirect, i)"%(planids[i], planids[i]))
        if no_beams == True:
            continue
        print("Setting plan file name:")
        
        #exec("tempmetainstuid = plands_%s.file_meta.MediaStorageSOPInstanceUID"%planids[i])

        tempmetainstuid = plansopinstuid + "." + str(i)

        planfilename = 'RP.' + tempmetainstuid + '.dcm'

        print("Plan file name: " + planfilename)
        planfilepath=Outputf + patientfolder + "/" + planfilename

        print(planfilepath)
        print("\n Saving plan file \n")
        exec("plands_%s.save_as(planfilepath)"%(planids[i]))
    #os.rename(Outputf+'%s'% patientfolder, Outputf+'%s,%s,%s'%(lname,fname,pid))
    #print("\n \n Current software versions found: \n")
    #for ver in softwarev:
        #print(ver)
    return softwarev
####################################################################################################################################################
####################################################################################################################################################



####################################################################################################################################################
# Function: initglobalvars
# This function simply resets all of the global variables to empty values
####################################################################################################################################################
def initglobalvars():
    global ROI_COUNT
    global structfilename
    global SeriesUID
    global StudyInstanceUID
    global FrameUID
    global ClassUID
    global patientname
    global dob
    global pid
    global imageslice
    global imageuid
    global Colors  
    global patient_sex
    global study_date
    global study_time
    global model
    global physician
    global sid
    global isocenter
    global ctcenter
    global descrip
    global plancount
    global plannamelist
    global planids
    global randval
    global currentdate
    global currenttime
    global doserefpt
    global patient_position
    global xshift
    global yshift
    global zshift
    global lname
    global fname
    global patientfolder
    global structsopinstuid
    global structseriesinstuid
    global plansopinstuid
    global planseriesinstuid
    global doseseriesuid
    global doseinstuid
    global planfilename
    global dosexdim
    global doseydim
    global dosezdim
    global doseoriginx
    global doseoriginy
    global doseoriginz
    global beamdosefiles
    global pixspacingx
    global pixspacingy
    global pixspacingz
    global posrefind
    global image_orientation
    global imagesetnumber
    global point_names
    global point_values
    global numfracs
    global flag_nobinaryfile
    global flag_noimages
    global no_setup_file
    global no_beams
    global softwarev
    global slicethick
    global x_dim
    global y_dim
    global z_dim
    global xpixdim
    global ypixdim

    ROI_COUNT = 0 #This value will represent the ROI that I'm currently looking at in file, will be incremented for each roi
    SeriesUID = "NA"
    StudyInstanceUID = "NA"
    FrameUID = "NA"
    ClassUID = "NA"
    patientname = ""
    dob = ""
    pid = ""
    imageslice = []
    imageuid = []
    Colors = [['255','0','0'],['255','20','147'],['0','0','255'],['0','255','0'],['125','38','205'],['255','255','0'],['255','140','0'],['0','100','0'],['0','191','255'],['255','192','203'],['72','209','204'],['139','69','19'],['255','193','37'],['221','160','221'],['107','142','35'],['142','35','35'],['245','204','176'],['191','239','255'],['139','28','98'],['255','99','71']] #red, pink, blue, green, purple, yellow, orange, dark green, sky blue, light pink, Turquois, brown, gold,lightpurple, olive, brick, peach?, light blue, maroon, tomato  
    patient_sex = ""
    study_date = ""
    study_time = ""
    model = ""
    physician = ""
    sid = ""
    isocenter = []
    ctcenter = []
    descrip = ""
    plancount = 0
    plannamelist = []
    planids = []
    randval = randint(0,999)
    currentdate = time.strftime("%Y%m%d")
    currenttime = time.strftime("%H%M%S")
    doserefpt = []
    patient_position = ""
    xshift = 0
    yshift = 0
    zshift = 0
    lname = ""
    fname = ""
    patientfolder = ""
    structsopinstuid = ''
    structseriesinstuid = ''
    plansopinstuid = ''
    planseriesinstuid = ''
    doseseriesuid = ''
    doseinstuid = ''
    planfilename = ''
    dosexdim = 0
    doseydim = 0
    dosezdim = 0
    doseoriginx = ""
    doseoriginy = ""
    doseoriginz = ""
    beamdosefiles = []
    pixspacingx = ""
    pixspacingy = ""
    pixspacingz = ""
    posrefind = ""
    image_orientation = []
    imagesetnumber = ""
    point_names = []
    point_values = []
    numfracs = ""
    flag_nobinaryfile = False
    flag_noimages = False
    no_setup_file = False
    no_beams = False
    softwarev = ""
    slicethick = 0
    x_dim = 0
    y_dim = 0
    z_dim = 0
    xpixdim = 0
    ypixdim = 0
####################################################################################################################################################
####################################################################################################################################################

####################################################################################################################################################
#    function: convertimages
#    The purpose of this function is to read in the image DICOM files
#    and to change the patients name to match the name of the patient
#    in the pinnacle files, also fills the list values for slicelocation
#    and UID. This function needs to be run, even if image files already converted
####################################################################################################################################################
def convertimages():
    print("Converting image patient name, birthdate and id to match pinnacle\n")
    global patientname
    global pid
    global dob
    global FrameUID
    global imageslice
    global SeriesUID
    global StudyInstanceUID
    global imageuid
    global patientfolder
    global posrefind
    global imagesetnumber
    global image_orientation
    global flag_noimages

    if not os.path.exists("%s%s/ImageSet_%s.DICOM"%(Inputf,patientfolder, imagesetnumber)):
        #Image set folder not found, need to ignore patient
        #Will want to call a function to be written that will create image set files from the condensed pixel data file
        print("Image files do not exist. Creating image files")
        createimagefiles()
        return
    for file in os.listdir("%s%s/ImageSet_%s.DICOM"%(Inputf,patientfolder, imagesetnumber)):
        if file == '11026.1.img':
            continue
        imageds = pydicom.read_file("%s%s/ImageSet_%s.DICOM/%s"%(Inputf, patientfolder, imagesetnumber, file), force=True)
        imageds.PatientsName = patientname
        imageds.PatientID = pid
        imageds.PatientsBirthDate = dob
        imageslice.append(imageds.SliceLocation)
        imageuid.append(imageds.SOPInstanceUID)
        image_orientation = imageds.ImageOrientationPatient
        tempinstuid = imageds.SOPInstanceUID
        posrefind = imageds.PositionReferenceIndicator
        imageds.SOPInstanceUID = tempinstuid
        imageds.FrameOfReferenceUID = FrameUID
        imageds.StudyInstanceUID = StudyInstanceUID
        imageds.SeriesInstanceUID = SeriesUID
        file_meta = Dataset()
        file_meta.TransferSyntaxUID = GTransferSyntaxUID
        file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.2'
        file_meta.MediaStorageSOPInstanceUID = tempinstuid
        file_meta.ImplementationClassUID = gImplementationClassUID
        imageds.file_meta = file_meta
        preamble = getattr(imageds, "preamble", None)
        if not preamble:
            preamble = b'\x00'*128
        currfile = DicomFile(Outputf+"%s/CT.%s.dcm"%(patientfolder, tempinstuid), 'wb')
        currfile.write(preamble)
        currfile.write(b'DICM')
        pydicom.write_file(Outputf+"%s/CT.%s.dcm"%(patientfolder,tempinstuid), imageds, False)
        print("Current image: ", file)
####################################################################################################################################################
####################################################################################################################################################                      
      

####################################################################################################################################################
# Function: createimagefiles()
# This function will create dicom image files for each slice using the condensed pixel data from file ImageSet_%s.img
####################################################################################################################################################
def createimagefiles():
    global slicethick
    global x_dim
    global y_dim
    global z_dim
    global xpixdim
    global ypixdim
    global patientname
    global pid
    global dob
    global FrameUID
    global imageslice
    global SeriesUID
    global StudyInstanceUID
    global imageuid
    global patientfolder
    global posrefind
    global imagesetnumber
    global image_orientation
    
    currentpatientposition = getheaderinfo()
    if os.path.isfile("%s%s/ImageSet_%s.img"%(Inputf, patientfolder, imagesetnumber)):
        allframeslist = []
        pixel_array = np.fromfile("%s%s/ImageSet_%s.img"%(Inputf, patientfolder, imagesetnumber), dtype = np.short)
        for i in range(0, int(z_dim)): # will loop over every frame
            frame_array = pixel_array[i*int(x_dim)*int(y_dim):(i+1)*int(x_dim)*int(y_dim)]
            allframeslist.append(frame_array)
            """frame_array = np.array([])
            temp_frame_array = pixel_array[i*int(x_dim)*int(y_dim):(i+1)*int(x_dim)*int(y_dim)]
            for j in range(0, int(y_dim)):
                temprow = temp_frame_array[j*int(x_dim):(j+1)*int(x_dim)][::-1]
                frame_array = np.append(frame_array, temprow)
            allframeslist.append(frame_array)
"""
    print("Length of frames list: " + str(len(allframeslist)))
    with open("%s%s/ImageSet_%s.ImageInfo"%(Inputf, patientfolder, imagesetnumber), 'rt', encoding='latin1') as f:
        image_info = f.readlines()
        curframe = 0
        for i, line in enumerate(image_info, 0):
            if "ImageInfo ={" in line:
                sliceloc = -float(re.findall(r"[-+]?\d*\.\d+|\d+", image_info[i + 1])[0])*10
                instuid = re.findall(r'"([^"]*)"', image_info[i + 8])[0]
                seriesuid = re.findall(r'"([^"]*)"', image_info[i + 4])[0]
                classuid = re.findall(r'"([^"]*)"', image_info[i + 7])[0]
                frameuid = re.findall(r'"([^"]*)"', image_info[i + 6])[0]
                studyinstuid = re.findall(r'"([^"]*)"', image_info[i + 5])[0]
                slicenum = int(re.findall(r"[-+]?\d*\.\d+|\d+", image_info[i + 3])[0])
                dateofscan, timeofscan = getdateandtime()
                
                file_meta = Dataset()
                file_meta.MediaStorageSOPClassUID = classuid
                file_meta.MediaStorageSOPInstanceUID = instuid
                file_meta.ImplementationClassUID = gImplementationClassUID #this value remains static since implementation for creating file is the same
                ds = FileDataset(planfilename, {}, file_meta=file_meta, preamble=b'\x00'*128)

                ds.SpecificCharacterSet = "ISO_IR 100"
                ds.ImageType = ['ORIGINAL', 'PRIMARY', 'AXIAL']
                ds.AccessionNumber = ''
                ds.SOPClassUID = classuid
                ds.SOPInstanceUID = instuid
                ds.StudyDate = dateofscan
                ds.SeriesDate = dateofscan
                ds.AcquisitionDate = dateofscan
                ds.ContentDate = dateofscan
                ds.AcquisitionTime = timeofscan
                ds.Modality = "CT" # Also should come from header file, but not always present
                ds.Manufacturer = "GE MEDICAL SYSTEMS" #This should come from Manufacturer in header, but for some patients it isn't set?? 
                ds.StationName = "CT"
                ds.PatientsName = patientname
                ds.PatientID = pid
                ds.PatientsBirthDate = dob
                ds.BitsAllocated = 16
                ds.BitsStored = 16
                ds.HighBit = 15
                ds.PixelRepresentation = 1
                ds.RescaleIntercept = -1024
                #ds.RescaleIntercept = 0.0
                ds.RescaleSlope = 1.0
                # ds.kvp = ?? This should be peak kilovoltage output of x ray generator used
                ds.PatientPosition = currentpatientposition
                ds.DataCollectionDiameter = xpixdim*float(x_dim)  # this is probably x_pixdim * xdim = y_pixdim * ydim
                ds.SpatialResolution = 0.35#???????
                #ds.DistanceSourceToDetector = #???
                #ds.DistanceSourceToPatient = #????
                ds.GantryDetectorTilt = 0.0 #??
                ds.TableHeight = -158.0#??
                ds.RotationDirection = "CW"#???
                ds.ExposureTime = 1000#??
                ds.XRayTubeCurrent = 398#??
                ds.GeneratorPower = 48#??
                ds.FocalSpots = 1.2#??
                ds.ConvolutionKernel = "STND" #????
                ds.SliceThickness = slicethick
                ds.NumberOfSlices =   int(z_dim)
                #ds.StudyInstanceUID = studyinstuid
                #ds.SeriesInstanceUID = seriesuid
                ds.FrameOfReferenceUID = FrameUID
                ds.StudyInstanceUID = StudyInstanceUID
                ds.SeriesInstanceUID = SeriesUID
                ds.InstanceNumber = slicenum # problem, some of these are repeated in image file so not sure what to do with that
                ds.ImagePositionPatient = [-xpixdim*float(x_dim)/2, -ypixdim*float(y_dim)/2, sliceloc]
                if "HFS" in currentpatientposition or "FFS" in currentpatientposition:
                    ds.ImageOrientationPatient = [1.0, 0.0, 0.0, 0.0, 1.0, -0.0] 
                elif "HFP" in currentpatientposition or "FFP" in currentpatientposition:
                    ds.ImageOrientationPatient = [-1.0, 0.0, 0.0, 0.0, -1.0, -0.0]
                ds.PositionReferenceIndicator = "LM" #???
                ds.SliceLocation = sliceloc
                ds.SamplesPerPixel = 1
                ds.PhotometricInterpretation = "MONOCHROME2"
                ds.Rows = int(x_dim)
                ds.Columns = int(y_dim)
                ds.PixelSpacing = [xpixdim, ypixdim]

                

                #ds.PixelData = allframeslist[curframe]
                #ds.PixelData = allframeslist[slicenum - 1]
                ds.PixelData = allframeslist[curframe].tostring()

                imageslice.append(sliceloc)
                imageuid.append(instuid)
                image_orientation = ds.ImageOrientationPatient
                posrefind = ds.PositionReferenceIndicator
                print("Creating image: " + Outputf + "%s/CT.%s.dcm"%(patientfolder, instuid))
                ds.save_as(Outputf + "%s/CT.%s.dcm"%(patientfolder, instuid))
                curframe = curframe + 1
####################################################################################################################################################
####################################################################################################################################################


####################################################################################################################################################
# Function: getheaderinfo
# This function will only be called in cases where image files do not already exist
####################################################################################################################################################
def getheaderinfo():
    global slicethick
    global x_dim
    global y_dim
    global z_dim
    global xpixdim
    global ypixdim
    temp_pos = ""
    with open("%s%s/ImageSet_%s.header"%(Inputf, patientfolder, imagesetnumber), "rt", encoding='latin1') as f2:
        for line in f2:
            #print("line in header: " + line)
            if "x_dim =" in line:
                x_dim = (line.split(" ")[-1]).replace(';','').replace('\n', '')
            if "y_dim =" in line:
                y_dim = (line.split(" ")[-1]).replace(';','').replace('\n', '')
            if "x_pixdim =" in line:
                xpixdim = float((line.split(" ")[-1]).replace(';',''))*10
            if "y_pixdim =" in line:
                ypixdim = float((line.split(" ")[-1]).replace(';',''))*10
            if "x_start =" in line and "index" not in line:
                xstart = float((line.split(" ")[-1]).replace(';',''))
                print("xstart = ", xstart)
            if "y_start =" in line:
                ystart = float((line.split(" ")[-1]).replace(';',''))
            if "z_dim =" in line:
                z_dim = (line.split(" ")[-1]).replace(';','').replace('\n', '')
            if "z_pixdim =" in line:
                slicethick = float((line.split(" ")[-1]).replace(';',''))*10
            if "z_start =" in line and "index" not in line:
                zstart = float((line.split(" ")[-1]).replace(';',''))
            if "patient_position" in line:
                temp_pos = (line.split(" ")[-1]).replace("\n","")
                print("Patient_position is: " + temp_pos)
    return temp_pos
####################################################################################################################################################
####################################################################################################################################################




####################################################################################################################################################
# Function: getdateandtime
# Will read ImageSet_%s.ImageSet file to get date and time of CT image aquisition, only used in cases where image files have not been created
####################################################################################################################################################
def getdateandtime():
    with open("//Testfile", "rt", encoding='latin1') as g:
        for line in g:
            if "ScanTimeFromScanner" in line:
                dateandtimestring = re.findall(r'"([^"]*)"', line)[0]
                dateandtimelist = dateandtimestring.split(' ')
                date = dateandtimelist[0].replace("-", "")
                time = dateandtimelist[1].replace(":","")
                return date, time
####################################################################################################################################################
####################################################################################################################################################

####################################################################################################################################################
#    function: readImageInfo
#    Reads in the file ImageSet_0.ImageInfo to get uid values
#    Saves the general UIDs as global variables 
####################################################################################################################################################
def readImageInfo():
    print("Reading image information for all image files\n")
    global SeriesUID
    global StudyInstanceUID
    global FrameUID
    global ClassUID
    global patientfolder
    global randval
    global imagesetnumber
    #print("Path to image info file: " + "%s%s/ImageSet_%s.ImageInfo"%(Inputf, patientfolder, imagesetnumber))
    if not os.path.exists("%s%s/ImageSet_%s.ImageInfo"%(Inputf, patientfolder, imagesetnumber)):
        #print("Leaving readImageInfo before getting info")
        return
    with open("%s%s/ImageSet_%s.ImageInfo"%(Inputf, patientfolder, imagesetnumber), 'rt', encoding='latin1') as f1:
        for line in f1:
            #print("For loop in readImageInfo")
            if "SeriesUID" in line:
                SeriesUID = re.findall(r'"([^"]*)"', line)[0]
                #print("setting series uid: " + str(SeriesUID))
                #SeriesUID = SeriesUID + "." + "0" + str(randval)
            if "StudyInstanceUID" in line:
                StudyInstanceUID = re.findall(r'"([^"]*)"', line)[0]
                #print("setting study uid: " + str(StudyInstanceUID))
                #StudyInstanceUID = StudyInstanceUID + "." + "0" + str(randval)
            if "FrameUID" in line:
                FrameUID = re.findall(r'"([^"]*)"', line)[0]
                #print("setting frame uid: " + str(FrameUID))
               # FrameUID = FrameUID[:-4] + "." + "0" + str(randval)
            if "ClassUID" in line:
                ClassUID = re.findall(r'"([^"]*)"', line)[0]
                #print("setting class uid: " + str(ClassUID))
                #ClassUID = ClassUID + "." + "0" + str(randval)
####################################################################################################################################################
####################################################################################################################################################


####################################################################################################################################################
# Creating a data structure to write to rt struct dicom file
# Based off example file write_new.py from C:\Python27\Lib\site-packages\dicom\examples\write_new.py
# returns data structure ds
####################################################################################################################################################
def createstructds():
    print("Creating Data structure")
    global structfilename
    global structsopinstuid
    # Populate required values for file meta information
    file_meta = Dataset()
    file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.481.3' #RT Structure Set Storage
    file_meta.MediaStorageSOPInstanceUID = structsopinstuid
    structfilename="RS."+structsopinstuid+".dcm"
    file_meta.ImplementationClassUID = gImplementationClassUID #this value remains static since implementation for creating file is the same
    # Create the FileDataset instance (initially no data elements, but file_meta supplied)
    ds = FileDataset(structfilename, {}, file_meta=file_meta, preamble=b'\x00'*128)
    #print(file_meta.preamble)
    return ds
####################################################################################################################################################
####################################################################################################################################################



####################################################################################################################################################
# Function: getstructshift()
# Purpose: reads in values from ImageSet_0.header to get x and y shift
####################################################################################################################################################
def getstructshift():
    global xshift
    global yshift
    global zshift
    global patient_position
    global imagesetnumber
    with open("%s%s/ImageSet_%s.header"%(Inputf, patientfolder, imagesetnumber), "rt", encoding='latin1') as f2:
        for line in f2:
            if "x_dim =" in line:
                x_dim = float((line.split(" ")[-1]).replace(';',''))
            if "y_dim =" in line:
                y_dim = float((line.split(" ")[-1]).replace(';',''))
            if "x_pixdim =" in line:
                xpixdim = float((line.split(" ")[-1]).replace(';',''))
            if "y_pixdim =" in line:
                ypixdim = float((line.split(" ")[-1]).replace(';',''))
            if "x_start =" in line and "index" not in line:
                xstart = float((line.split(" ")[-1]).replace(';',''))
                print("xstart = ", xstart)
            if "y_start =" in line:
                ystart = float((line.split(" ")[-1]).replace(';',''))
            if "z_dim =" in line:
                z_dim = float((line.split(" ")[-1]).replace(';',''))
            if "z_pixdim =" in line:
                zpixdim = float((line.split(" ")[-1]).replace(';',''))
            if "z_start =" in line and "index" not in line:
                zstart = float((line.split(" ")[-1]).replace(';',''))
    if patient_position == 'HFS':
        xshift = ((x_dim*xpixdim/2)+xstart)*10
        print("X shift = ", xshift)
        yshift = -((y_dim*ypixdim/2)+ystart)*10
        print("Y shift = ", yshift)
        zshift = -((z_dim*zpixdim/2)+zstart)*10
    elif patient_position == 'HFP':
        xshift = -((x_dim*xpixdim/2)+xstart)*10

        print("X shift = ", xshift)


        yshift = ((y_dim*ypixdim/2)+ystart)*10

        print("Y shift = ", yshift)
        zshift = -((z_dim*zpixdim/2)+zstart)*10
    elif patient_position == 'FFP':
        xshift = ((x_dim*xpixdim/2)+xstart)*10
        print("X shift = ", xshift)
        yshift = ((y_dim*ypixdim/2)+ystart)*10
        print("Y shift = ", yshift)
        zshift = ((z_dim*zpixdim/2)+zstart)*10
    elif patient_position == 'FFS':
        xshift = -((x_dim*xpixdim/2)+xstart)*10
        print("X shift = ", xshift)
        yshift = -((y_dim*ypixdim/2)+ystart)*10
        print("Y shift = ", yshift)
        zshift = ((z_dim*zpixdim/2)+zstart)*10

####################################################################################################################################################
####################################################################################################################################################


####################################################################################################################################################
# Function: createplands()
# creates a data structure for the rt plan file
# similar to createstructds but different UIDs
####################################################################################################################################################
def createplands(plannumber):
    print("Creating plan Data Structure\n")
    global planfilename
    global plansopinstuid
    file_meta = Dataset()
    file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.481.5' #RT Plan Storage
    file_meta.MediaStorageSOPInstanceUID = plansopinstuid + "." + str(plannumber)
    file_meta.ImplementationClassUID = gImplementationClassUID #this value remains static since implementation for creating file is the same
    ds = FileDataset(planfilename, {}, file_meta=file_meta, preamble=b'\x00'*128)
    return ds
####################################################################################################################################################
####################################################################################################################################################

####################################################################################################################################################
# Function to initialize data structure, sets Specific Character set,
# instance creation time and date, SOP Class UID, and SOP Instance UID
# also sets modality and data structure
####################################################################################################################################################
def initds(ds):
    print("initializing data structure\n")
    global imageuids    
    global SeriesUID
    global StudyInstanceUID
    global FrameUID
    global ClassUID
    global structsopinstuid
    global structseriesinstuid
    ds.SpecificCharacterSet = 'ISO_IR 100' # not sure what I want here, going off of template dicom file
    ds.InstanceCreationDate = time.strftime("%Y%m%d")
    ds.InstanceCreationTime = time.strftime("%H%M%S")
    ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.481.3'
    ds.SOPInstanceUID = structsopinstuid
    ds.Modality = 'RTSTRUCT'
    ds.AccessionNumber = ""
    ds.Manufacturer = Manufacturer #from sample dicom file, maybe should change?
    ds.StationName = "adacp3u7" # not sure where to get information for this element can find this and read in from 
    #ds.ManufacturersModelName = 'Pinnacle3'
    ReferencedStudy1 = Dataset()
    ds.ReferencedStudySequence.append(ReferencedStudy1)
    ds.ReferencedStudySequence[0].ReferencedSOPClassUID = '1.2.840.10008.3.1.2.3.2' #Study Component Management SOP Class (chosen from template)
    ds.ReferencedStudySequence[0].ReferencedSOPInstanceUID = StudyInstanceUID
    ds.StudyInstanceUID = StudyInstanceUID
    print("Setting structure file study instance: " + str(StudyInstanceUID))
    ds.SeriesInstanceUID = structseriesinstuid
    return ds
####################################################################################################################################################
####################################################################################################################################################


####################################################################################################################################################
# function to read in patient info from Patient text file
####################################################################################################################################################
def readpatientinfo(ds):
    print ("Reading patient information\n")
    mname = ""
    flag_first = True
    flag_time = True
    flag_stid = True
    global imageuids    
    global SeriesUID
    global StudyInstanceUID
    global FrameUID
    global ClassUID
    global patientname
    global dob
    global pid
    global patient_sex
    global study_time
    global study_date
    global model
    global physician
    global sid
    global descrip
    global plancount
    global plannamelist
    global patientfolder
    global lname
    global fname
    global imagesetnumber
    global softwarev
    #global listofversions
    with open("%s%s/Patient"%(Inputf, patientfolder), "rt", encoding='latin1') as g: 
        for line in g:
            if "PatientID =" in line:
                pid = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
                ds.PatientID = pid #may want to change pid to be value of medical record number.
            if "LastName = " in line:
                lname = re.findall(r'"([^"]*)"', line)[0]
                lname = lname.replace(' (restored)', '')
                lname = lname.replace('\\', '')
                lname = lname.replace('/', '')
            if "FirstName =" in line:
                fname = re.findall(r'"([^"]*)"', line)[0]
                fname = fname.replace(" ", "")
                fname = fname.replace("\\", '')
                fname = fname.replace("/", '')
            if "MiddleName = " in line:
                mname = re.findall(r'"([^"]*)"', line)[0]
                mname = mname.replace("\\", '')
                mname = mname.replace('/', '')
            if "MedicalRecordNumber =" in line:
                medrecnum = re.findall(r'"([^"]*)"', line)[0] 
            if "ReferringPhysician = " in line:
                refphys = re.findall(r'"([^"]*)"', line)[0]
                ds.ReferringPhysiciansName = refphys
            if "RadiationOncologist = " in line:
                physician = re.findall(r'"([^"]*)"', line)[0]
                ds.PhysiciansOfRecord = physician
            if "Comment = " in line and flag_first:
                descrip = re.findall(r'"([^"]*)"', line)[0]
                ds.StudyDescription = descrip
            if "Gender = " in line:
                gen = re.findall(r'"([^"]*)"', line)[0]
                if "Male" in gen:
                    patient_sex = 'M'
                elif "Female" in gen:
                    patient_sex = 'F'
                ds.PatientsSex = patient_sex
            if "DateOfBirth =" in line:
                dobstr = re.findall(r'"([^"]*)"', line) #gets birthday string with numbers and dashes
                if '-' in dobstr[0]:
                    dob_list = dobstr[0].split('-')
                elif '/' in dobstr[0]:
                    dob_list = dobstr[0].split('/')
                else:
                    dob_list = dobstr[0].split(' ')
                dob = ""
                for num in dob_list:
                    if len(num) == 1:
                        num = '0' + num
                    if num == dob_list[-1] and len(num) == 2:
                        num = "19" + num
                    dob = dob + num
                ds.PatientsBirthDate = dob
            if "ImageSetList ={" in line:
                flag_first = False
            if "PlanName = " in line:
                plancount = plancount + 1
                plannamelist.append(re.findall(r'"([^"]*)"', line)[0])
                ds.StructureSetLabel = plannamelist[plancount - 1]
            if "PrimaryCTImageSetID =" in line:
                imagesetnumber = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
                print("Image set number: " + imagesetnumber)
            if "    PlanID =" in line:
                planids.append(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])
            if "    StudyID = " in line and flag_stid:
                sid = re.findall(r'"([^"]*)"', line)[0]
                print ("Study id: ", sid)
                ds.StudyID = sid
                flag_stid = False
            if "WriteTimeStamp = " in line and flag_time:
                dateandtime = re.findall(r'"([^"]*)"', line)[0]
                arr = dateandtime.split()
                date = arr[0].replace('-', '')
                time = arr[1].replace(':', '')
                ds.StructureSetDate = date
                ds.StructureSetTime = time
                study_date = date
                study_time = time
                ds.StudyDate = date
                ds.StudyTime = time
                flag_time = False
            if "ToolType =" in line:
                model = re.findall(r'"([^"]*)"', line)[0]
                ds.ManufacturersModelName = model
            if "PinnacleVersionDescription" in line:
                softwarev = re.findall(r'"([^"]*)"', line)[0]
                ds.SoftwareVersions = softwarev
                #if listofversions == []:
                   # listofversions.append(softwarev)
                #flagsame = False
                #for ver in listofversions:
                    #if softwarev == ver:
                        #flagsame = True
                #f flagsame == False:
                    #listofversions.append(softwarev)
    ds.StructureSetName = 'POIandROI'
    ds.SeriesNumber = '1'
    patientname = lname + "^" + fname + "^" + mname + "^"
    ds.PatientsName = patientname
    return ds
####################################################################################################################################################
####################################################################################################################################################


####################################################################################################################################################
#  Function that read in the first ROIs (reference points) from plan.Points
# Takes in data structure
####################################################################################################################################################
def readpoints(ds, planfolder):
    print("Reading in the points\n")
    global ROI_COUNT 
    global SeriesUID
    global StudyInstanceUID
    global FrameUID
    global ClassUID
    global imageslice
    global imageuid
    global Colors
    global isocenter
    global patientfolder
    global doserefpt
    global ctcenter
    global xshift
    global yshift
    global patient_position
    global point_values
    global point_names
    global FrameUID
    with open("%s%s/%s/plan.Points"%(Inputf, patientfolder, planfolder), "rt", encoding='latin1') as e:
        for num, line in enumerate(e,1):
            if "  Name = " in line:
                ROI_COUNT = ROI_COUNT + 1
                exec("ROIContour%d = Dataset()"%(ROI_COUNT))
                exec("ds.ROIContourSequence.append(ROIContour%d)"%ROI_COUNT)
                exec("ROIContour%d.ReferencedROINumber = str(%d)"%(ROI_COUNT, ROI_COUNT))
                exec("StructureSetROI%d = Dataset()"%(ROI_COUNT))
                exec("ds.StructureSetROISequence.append(StructureSetROI%d)"%(ROI_COUNT))
                exec("RTROIObservations%d = Dataset()"%(ROI_COUNT))
                exec("ds.RTROIObservationsSequence.append(RTROIObservations%d)"%ROI_COUNT)
                ds.StructureSetROISequence[ROI_COUNT - 1].ROINumber = ROI_COUNT
                refptname = re.findall(r'"([^"]*)"', line)[0]    
                ds.StructureSetROISequence[ROI_COUNT - 1].ROIName = refptname
                point_names.append(refptname)
                ds.ROIContourSequence[ROI_COUNT - 1].ContourSequence = Sequence()
                ds.StructureSetROISequence[ROI_COUNT - 1].ROIGenerationAlgorithm = 'SEMIAUTOMATIC' #Not sure what this is for, just basing off template, should look into further
                ds.StructureSetROISequence[ROI_COUNT - 1].ReferencedFrameofReferenceUID = FrameUID
                ds.ROIContourSequence[ROI_COUNT -1].ROIDisplayColor = Colors[0]
                refpoint = []
            if "XCoord =" in line:
                x = (line.split(" ")[-1]).replace(';','')
                if patient_position == 'HFS' or patient_position == 'FFP':
                    x = str(float(x)*10)
                elif patient_position == 'HFP' or patient_position == 'FFS':
                    x = str(-float(x)*10)
                refpoint.append(x)
            if "YCoord =" in line:
                y = (line.split(" ")[-1]).replace(';','')
                if patient_position == 'HFS' or patient_position == 'FFS':
                    y = str(-float(y)*10)
                elif patient_position == 'HFP' or patient_position == 'FFP':
                    y = str(float(y)*10)
                refpoint.append(y)
            if "ZCoord =" in line:
                z = (line.split(" ")[-1]).replace(';','')
                #print("The line below is the z value:")
                #print(z)
                #print("The line above is the z value")
                if patient_position == 'HFS' or patient_position == 'HFP':
                    z = str(-float(z)*10)
                elif patient_position == 'FFS' or patient_position == 'FFP':
                    z = str(float(z)*10)
                refpoint.append(z)
            if "LastModifiedTimeStamp" in line: #this is the last line for the points
                ds.RTROIObservationsSequence[ROI_COUNT -1].ObservationNumber = ROI_COUNT
                ds.RTROIObservationsSequence[ROI_COUNT -1].ReferencedROINumber = ROI_COUNT
                ds.RTROIObservationsSequence[ROI_COUNT -1].RTROIInterpretedType = 'MARKER'
                ds.RTROIObservationsSequence[ROI_COUNT -1].ROIInterpreter = ""
                Contour1 = Dataset()
                ds.ROIContourSequence[ROI_COUNT - 1].ContourSequence.append(Contour1)
                ds.ROIContourSequence[ROI_COUNT - 1].ContourSequence[0].ContourData = refpoint
                point_values.append(refpoint)
                print("refpoint:", refpoint)
                #if isocenter == []:
                #   isocenter = refpoint
                if "Iso" in refptname or "isocenter" in refptname or "isocentre" in refptname:
                    isocenter = refpoint
                if "CT Center" in refptname or "ct center" in refptname or "ct centre" in refptname:
                    ctcenter = refpoint
                if "drp" in refptname or "DRP" in refptname:
                    doserefpt = refpoint
                print("isocenter:", isocenter)
                ds.ROIContourSequence[ROI_COUNT - 1].ContourSequence[0].ContourGeometricType = 'POINT'
                ds.ROIContourSequence[ROI_COUNT - 1].ContourSequence[0].NumberofContourPoints = 1
                ds.ROIContourSequence[ROI_COUNT - 1].ContourSequence[0].ContourImageSequence = Sequence()
                ContourImage1 = Dataset()
                ds.ROIContourSequence[ROI_COUNT - 1].ContourSequence[0].ContourImageSequence.append(ContourImage1)
                closestvalue = abs(float(imageslice[0]) - float(refpoint[-1]))
                closestlocation = 0
                match = False
                for i, s in enumerate(imageslice,0):
                    #print("finding corresponding image\n")
                    if abs(float(s) - (float(refpoint[-1]))) < 0.01: #making this the tolerance
                        ds.ROIContourSequence[ROI_COUNT - 1].ContourSequence[0].ContourImageSequence[0].ReferencedSOPClassUID = '1.2.840.10008.5.1.4.1.1.2'
                        ds.ROIContourSequence[ROI_COUNT - 1].ContourSequence[0].ContourImageSequence[0].ReferencedSOPInstanceUID = imageuid[i]
                        match = True
                    else:
                        if abs(float(s) - (float(refpoint[-1]))) < closestvalue:
                            closestvalue = abs(float(s) - (float(refpoint[-1])))
                            closestlocation = i
                if not match:
                    ds.ROIContourSequence[ROI_COUNT - 1].ContourSequence[0].ContourImageSequence[0].ReferencedSOPClassUID = '1.2.840.10008.5.1.4.1.1.2'
                    ds.ROIContourSequence[ROI_COUNT - 1].ContourSequence[0].ContourImageSequence[0].ReferencedSOPInstanceUID = imageuid[closestlocation]

    if len(isocenter) < 2:
        isocenter = ctcenter
        print("Isocenter not located, setting to ct center: ", str(isocenter))
    if len(isocenter) < 2:
        #print("Isocenter still not located, setting to point with center in name, if not, with iso in name")
        temp_point1 = []
        temp_point2 = []
        for j, pointnames in enumerate(point_names):
            if "center" in pointnames:
                temp_point1 = point_values[j]
                #print("setting to: " + str(temp_point1))
            elif "iso" in pointnames:
                temp_point2 = point_values[j]
                #print("setting to: " + str(temp_point2))
        if len(temp_point1) > 1:
            isocenter = temp_point1
            #print("setting iso: " + str(isocenter))
        elif len(temp_point2) > 1:
            isocenter = temp_point2
            #print("setting iso: " + str(isocenter))
        else:
            if  len(ds.ROIContourSequence) > 0: 
                isocenter = point_values[0] # setting to first point if isocenter or ct center not found
                #print("setting iso to actual value: " + str(isocenter))
    print("isocenter before loop to apply shifts to contour sequence points: " + str(isocenter))
    for enteredpoints in ds.ROIContourSequence:
        #print("In loop applying shifts: isocenter:" + str(isocenter) )
        enteredpoints.ContourSequence[0].ContourData[0] = str(float(enteredpoints.ContourSequence[0].ContourData[0]) - xshift)
        enteredpoints.ContourSequence[0].ContourData[1] = str(float(enteredpoints.ContourSequence[0].ContourData[1]) - yshift)
        #enteredpoints.ContourSequence[0].ContourData[2] = str(float(enteredpoints.ContourSequence[0].ContourData[2]) - float(isocenter[2]))  
        #print("bottom of loop applying shifts isocenter:" + str(isocenter))  
    print("end of read points isocenter:" + str(isocenter))
    return ds            
####################################################################################################################################################
####################################################################################################################################################


####################################################################################################################################################
# Function to read in plan.roi and get contour information
# will take in data structure.
####################################################################################################################################################
def readroi(ds, planfolder):
    print("Reading in roi file\n")
    global ROI_COUNT   
    global SeriesUID
    global StudyInstanceUID
    global FrameUID
    global ClassUID
    global imageslice
    global imageuid
    global Colors
    global patientfolder
    global isocenter
    global xshift
    global yshift   
    global patient_position
    points = []
    flag_points = False # bool value to tell me if I want to read the line in as point values
    prevroi = ROI_COUNT
    with open("%s%s/%s/plan.roi"%(Inputf, patientfolder, planfolder), "rt", encoding='latin1') as f:
        for num, line in enumerate(f, 1):
            if "};  // End of points for curve" in line: # this will tell me not to read in point values
                #all points for current curve saved until now. Here is where I need to add them to dicom file
                ds.ROIContourSequence[ROI_COUNT - 1].ContourSequence[int(curvenum) - 1].ContourData = points
                ds.ROIContourSequence[ROI_COUNT - 1].ContourSequence[int(curvenum) - 1].ContourImageSequence = Sequence()
                ContourImage1 = Dataset()
                flag_match = False
                ds.ROIContourSequence[ROI_COUNT - 1].ContourSequence[int(curvenum) - 1].ContourImageSequence.append(ContourImage1)
                closestvalue = abs(float(imageslice[0]) - (float(points[-1])))
                closestlocation = 0
                match = False
                for i, s in enumerate(imageslice,0):
                    #print("finding corresponding image\n")
                    if abs(float(s) - (float(points[-1]))) < 0.01: #making this the tolerance for a match
                        ds.ROIContourSequence[ROI_COUNT - 1].ContourSequence[int(curvenum) - 1].ContourImageSequence[0].ReferencedSOPClassUID = '1.2.840.10008.5.1.4.1.1.2'
                        ds.ROIContourSequence[ROI_COUNT - 1].ContourSequence[int(curvenum) - 1].ContourImageSequence[0].ReferencedSOPInstanceUID = imageuid[i]
                        match = True
                    else:
                        if abs(float(s) - (float(points[-1]))) < closestvalue:
                            closestvalue = abs(float(s) - (float(points[-1])))
                            closestlocation = i
                if not match:
                    ds.ROIContourSequence[ROI_COUNT - 1].ContourSequence[int(curvenum) - 1].ContourImageSequence[0].ReferencedSOPClassUID = '1.2.840.10008.5.1.4.1.1.2'
                    ds.ROIContourSequence[ROI_COUNT - 1].ContourSequence[int(curvenum) - 1].ContourImageSequence[0].ReferencedSOPInstanceUID = imageuid[closestlocation]
                del points[:]
                flag_points = False
            if flag_points:
                curr_points = line.split(' ')
                if patient_position == 'HFS':
                    curr_points = [str(float(curr_points[0])*10 - xshift), str(-float(curr_points[1])*10 - yshift), str(-float(curr_points[2])*10)]
                elif patient_position == 'HFP':
                    curr_points = [str(-float(curr_points[0])*10 - xshift), str(float(curr_points[1])*10 - yshift), str(-float(curr_points[2])*10)]
                elif patient_position == 'FFP':
                    curr_points = [str(float(curr_points[0])*10 - xshift), str(float(curr_points[1])*10 - yshift), str(float(curr_points[2])*10)]
                elif patient_position == 'FFS':
                    curr_points = [str(-float(curr_points[0])*10 - xshift), str(-float(curr_points[1])*10 - yshift), str(float(curr_points[2])*10)]
                points = points + curr_points
            if "Beginning of ROI" in line: # Start of ROI
                ROI_COUNT = ROI_COUNT + 1 #increment ROI_num because I've found a new ROI
                exec("ROIContour%d = Dataset()"%(ROI_COUNT))
                exec("ds.ROIContourSequence.append(ROIContour%d)"%ROI_COUNT)
                exec("ROIContour%d.ReferencedROINumber = str(%d)"%(ROI_COUNT, ROI_COUNT))
                exec("StructureSetROI%d = Dataset()"%(ROI_COUNT))
                exec("ds.StructureSetROISequence.append(StructureSetROI%d)"%(ROI_COUNT))
                exec("RTROIObservations%d = Dataset()"%(ROI_COUNT))
                exec("ds.RTROIObservationsSequence.append(RTROIObservations%d)"%ROI_COUNT)
                ds.StructureSetROISequence[ROI_COUNT - 1].ROINumber = ROI_COUNT
                ROIName = line[22:] # gets a string of ROI name
                ROIName = ROIName.replace('\n','')
                ds.StructureSetROISequence[ROI_COUNT - 1].ROIName = ROIName
                ds.StructureSetROISequence[ROI_COUNT - 1].ROIGenerationAlgorithm = 'SEMIAUTOMATIC'
                ds.StructureSetROISequence[ROI_COUNT - 1].ReferencedFrameofReferenceUID = FrameUID
                ds.ROIContourSequence[ROI_COUNT - 1].ContourSequence = Sequence()
                if ROI_COUNT - prevroi <= len(Colors):
                    ds.ROIContourSequence[ROI_COUNT -1].ROIDisplayColor = Colors[ROI_COUNT - prevroi - 1]
                    #print(ROI_COUNT - prevroi-1)
                else:
                    if ROI_COUNT - 1 - len(Colors) < len(Colors):
                        ds.ROIContourSequence[ROI_COUNT -1].ROIDisplayColor = Colors[ROI_COUNT - 1 - len(Colors)]
                        #print(ROI_COUNT - 1 - len(Colors))
                    elif ROI_COUNT - 1 - len(Colors) - len(Colors) < len(Colors):
                        ds.ROIContourSequence[ROI_COUNT -1].ROIDisplayColor = Colors[ROI_COUNT - 1 - len(Colors) - len(Colors)]
                        #print(ROI_COUNT - 1 - len(Colors) - len(Colors))
                    else:
                        ds.ROIContourSequence[ROI_COUNT -1].ROIDisplayColor = Colors[ROI_COUNT - 1 - len(Colors) - len(Colors) - len(Colors)]
                        #print(ROI_COUNT - 1 - len(Colors) - len(Colors) - len(Colors))
                print('ROI Number & Name: '+str(ROI_COUNT)+', '+ROIName)
            if "}; // End of ROI" in line: #end of ROI found
                #ROI_type = line[31:]
                #ROI_type = ROI_type.replace('\n','')
                ds.RTROIObservationsSequence[ROI_COUNT -1].ObservationNumber = ROI_COUNT
                ds.RTROIObservationsSequence[ROI_COUNT -1].ReferencedROINumber = ROI_COUNT
                if "PTV" in ROIName: 
                    ds.RTROIObservationsSequence[ROI_COUNT -1].RTROIInterpretedType = 'PTV'
                else:
                    ds.RTROIObservationsSequence[ROI_COUNT -1].RTROIInterpretedType = 'ORGAN'
                ds.RTROIObservationsSequence[ROI_COUNT -1].ROIInterpreter = ""
                #add to ROI observation sequence
            if "volume =" in line:
                vol = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
                ds.StructureSetROISequence[ROI_COUNT - 1].ROIVolume = vol
            if "//  Curve " in line: #found a curve
                curvenum = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
                exec("Contour%s = Dataset()"%curvenum)
                exec("ds.ROIContourSequence[%d - 1].ContourSequence.append(Contour%d)"%(ROI_COUNT, int(curvenum)))
            if "num_points =" in line:
                npts = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
                ds.ROIContourSequence[ROI_COUNT - 1].ContourSequence[int(curvenum) - 1].ContourGeometricType = 'CLOSED_PLANAR'
                ds.ROIContourSequence[ROI_COUNT - 1].ContourSequence[int(curvenum) - 1].NumberofContourPoints = npts
            if "points=" in line:
                flag_points = True
    return ds
####################################################################################################################################################
####################################################################################################################################################

####################################################################################################################################################
# Function: planinit()
# purpose: to fill in data basic elements for RT plan file
####################################################################################################################################################
def planinit(ds, planame, planandnum, plannumber):
    global patientname
    global dob
    global pid
    global patient_sex
    global plansopinstuid
    global study_time
    global study_date
    global StudyInstanceUID
    global model
    global physician
    global planseriesinstuid
    global sid
    global FrameUID
    global descrip
    global structsopinstuid
    global posrefind
    ds.SpecificCharacterSet = 'ISO_IR 100'
    ds.InstanceCreationDate = time.strftime("%Y%m%d")
    ds.InstanceCreationTime = time.strftime("%H%M%S")
    ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.481.5' #RT Plan Storage
    ds.SOPInstanceUID = plansopinstuid+"." + str(plannumber)
    ds.StudyDate = study_date
    ds.StudyTime = study_time
    ds.AccessionNumber = ''
    ds.Modality = 'RTPLAN'
    ds.Manufacturer =Manufacturer
    ds.OperatorsName = ""
    ds.ManufacturersModelName = model
    ds.SoftwareVersions = ['u\'9.0']
    ds.PhysiciansOfRecord = physician
    ds.PatientsName = patientname
    ds.PatientsBirthDate = dob
    ds.PatientID = pid
    ds.PatientsSex = patient_sex
    ds.StudyInstanceUID = StudyInstanceUID
    ds.SeriesInstanceUID = planseriesinstuid + "." + str(plannumber)
    ds.StudyID = sid
    ds.FrameofReferenceUID = FrameUID
    ds.PositionReferenceIndicator = posrefind
    ds.RTPlanLabel = planandnum + '.0' # may need to change this later
    ds.RTPlanName = planame
    ds.RTPlanDescription = descrip
    ds.RTPlanDate = study_date
    ds.RTPlanTime = study_time
    #ds.PlanIntent = "" #Not sure where to get this informationd, will likely be 'CURATIVE' or 'PALIATIVE'
    ds.RTPlanGeometry = 'PATIENT'
    #ds.DoseReferenceSequence = Sequence() #figure out what goes in DoseReferenceSequence... Should be like a target volume and reference point I think...
    #ds.ToleranceTableSequence = Sequence() #figure out where to get this information
    ds.FractionGroupSequence = Sequence()
    ds.BeamSequence = Sequence()
    ds.PatientSetupSequence = Sequence() #need one per beam
    ds.ReferencedStructureSetSequence = Sequence()
    ReferencedStructureSet1 = Dataset()
    ds.ReferencedStructureSetSequence.append(ReferencedStructureSet1)
    ds.ReferencedStructureSetSequence[0].ReferencedSOPClassUID = '1.2.840.10008.5.1.4.1.1.481.3'
    ds.ReferencedStructureSetSequence[0].ReferencedSOPInstanceUID = structsopinstuid
    ds.ApprovalStatus = 'UNAPPROVED' #find out where to get this information
    return ds
####################################################################################################################################################
####################################################################################################################################################


####################################################################################################################################################
# Function: getpatientsetup()
# purpose: Returns the value for patient setup
####################################################################################################################################################
def getpatientsetup(planfolder):
    global patientfolder
    global no_setup_file
    if not os.path.exists("%s%s/%s/plan.PatientSetup"%(Inputf, patientfolder, planfolder)):
        no_setup_file = True
        return
    with open("%s%s/%s/plan.PatientSetup"%(Inputf, patientfolder, planfolder), "rt", encoding='latin1') as f:
        for line in f:
            if "Position =" in line:
                pos = re.findall(r'"([^"]*)"', line)[0]
            if "Orientation =" in line:
                orient = re.findall(r'"([^"]*)"', line)[0]
        if "Head First" in orient:
            pat_pos = "HF"
        elif "Feet First" in orient:
            pat_pos = "FF"
        if "supine" in pos:
            pat_pos = pat_pos + "S"
        elif "prone" in pos:
            pat_pos = pat_pos + "P"
        elif "decubitus right" in pos or "Decuibitus Right" in pos:
            pat_pos = pat_pos + "DR"
        elif "decubitus left" in pos or "Decuibitus Left" in pos:
            pat_pos = pat_pos +"DL"
    return pat_pos
####################################################################################################################################################
####################################################################################################################################################


####################################################################################################################################################
# Function: readtrial()
# purpose: get Beam information from plan.Trial for RT plan file
####################################################################################################################################################
def readtrial(ds, planfolder, plannumber):
    #print("There is a problem somewhere in this function\n")
    global isocenter
    global xshift
    global yshift
    global zshift
    global patientfolder
    global doserefpt
    global patient_position
    global dosexdim
    global doseydim
    global dosezdim
    global doseoriginx
    global doseoriginy
    global doseoriginz
    global beamdosefiles
    global numfracs
    global pixspacingy
    global pixspacingx
    global pixspacingz
    global point_values
    global point_names
    global flag_nobinaryfile
    global no_beams
    global PDD15MV
    global PDD16MV
    global PDD10MV
    global PDD6MV
    print("Entering readtrial function, isocenter: " + str(isocenter))
    beamdoses = []
    beamdosefiles = []
    beamcount = 0
    MUlineflag = False
    nomachinename = True
    mlcleafpos = False
    flag_stepnshoot = False
    beginleafpoints = False
    currentmeterset = 0.0
    beginbeam = False
    ctrlptlist = False
    ctrlptmeterflag = False
    noname = True
    current_dosefile_num = ''
    beamenergies = []
    leafpositions1 = []
    leafpositions2 = []
    metersetweight = ['0']
    wedgeangles = []
    numctrlpts = 0
    totalleafpositions = [] #this list will have lists for all the control points
    currentcontrolpoint = 0
    countpoints = 0
    FractionGroup1 = Dataset() #I'm assuming here I only need one data set in fraction goup sequence
    ds.FractionGroupSequence.append(FractionGroup1)
    ds.FractionGroupSequence[0].ReferencedBeamSequence = Sequence()
    tempfile = open("%s%s/%s/plan.Trial"%(Inputf, patientfolder, planfolder), "rt", encoding='latin1')
    all_lines = tempfile.readlines() #this is a big waste of space, there is probably a better way to do this, but it will work for now
    tempfile.close()
    num_trials = all_lines.count("Trial ={\n")
    if num_trials > 1:
        #for i in range(0, num_trials):
            #if re.findall(r"[-+]?\d*\.\d+|\d+", next((s for s in all_lines if " UseTrialForTreatment" in s), None))[0] == '0':
                #linetostart = all_lines[1:].index("Trial ={\n") + 1
                #all_lines = all_lines[linetostart:]
        linetostart = all_lines[1:].index("Trial ={\n") + 1
        all_lines = all_lines[:linetostart] #take first trial for now
    #with open("%s%s/%s/plan.Trial"%(Inputf, patientfolder, planfolder), "rt", encoding='latin1') as h:
        #for linenum, line in enumerate(h,0):
    for linenum, line in enumerate(all_lines, 0):
        if "BeamList ={" in line and "};" in all_lines[linenum + 1]:
            #empty beam set, skip patient
            no_beams = True
            return
        if "DoseGrid .VoxelSize .X" in line:
            pixspacingx = float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])*10
        if "DoseGrid .VoxelSize .Y" in line:
            pixspacingy = float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])*10
        if "DoseGrid .VoxelSize .Z" in line:
            pixspacingz = float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])*10
        if "DoseGrid .Dimension .X" in line:
            dosexdim = int(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])
        if "DoseGrid .Dimension .Y" in line:
            doseydim = int(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])
        if "DoseGrid .Dimension .Z" in line:
            dosezdim = int(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])
        if "DoseGrid .Origin .X" in line:
            if patient_position == 'HFP' or patient_position == 'FFS':
                doseoriginx = str(-float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])*10 - xshift)
            elif patient_position == 'HFS' or patient_position == 'FFP':
                doseoriginx = str(float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])*10 - xshift)
             
        if "DoseGrid .Origin .Y" in line:
            if patient_position == 'HFS' or patient_position == 'FFS':
                doseoriginy = str(-float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])*10 - yshift)
            elif patient_position == 'HFP' or patient_position == 'FFP':
                doseoriginy = str(float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])*10 - yshift)
        if "DoseGrid .Origin .Z" in line:
            if patient_position == 'HFS' or patient_position == 'HFP':
                doseoriginz = str(-float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])*10)
            elif patient_position == 'FFS' or patient_position == 'FFP':
                doseoriginz = str(float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])*10)
        if "      NumberOfFractions =" in line:
            numfracs = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "      DoseVolume = " in line:
            current_dosefile_num = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
            if int(current_dosefile_num) < 10:
                current_dosefile_num = "00" + current_dosefile_num
            elif int(current_dosefile_num) < 100:
                current_dosefile_num = "0" + current_dosefile_num 
            beamdosefiles.append(current_dosefile_num)
            print('Reading file: plan.Trail.binary.'+str(current_dosefile_num))
            print('Number of dose files read: '+str(len(beamdosefiles)))
        if "Beam ={" in line and 'Proton' not in line:
            #print("Line that indicates beam information\n")
            #new beam
            MUlineflag = False
            nomachinename = True
            noname = True
            countpoints = 0
            currentcontrolpoint = 0
            numwedges = 0
            beginbeam = True
            wedgeflag = False
            beamcount = beamcount + 1
            del totalleafpositions
            totalleafpositions = []
            del leafpositions1
            del leafpositions2
            leafpositions1 = []
            leafpositions2 = []
            x1 = ""
            x2 = ""
            y1 = ""
            y2 = ""
            exec("ReferencedBeam%d = Dataset()"%beamcount)
            exec("ds.FractionGroupSequence[0].ReferencedBeamSequence.append(ReferencedBeam%d)"%beamcount)
             
            ds.FractionGroupSequence[0].ReferencedBeamSequence[beamcount - 1].ReferencedBeamNumber = beamcount
            exec("Beam%d = Dataset()"%beamcount)
            exec("ds.BeamSequence.append(Beam%d)"%beamcount)
            ds.BeamSequence[beamcount - 1].Manufacturer = Manufacturer #figure out what to put here
            ds.BeamSequence[beamcount - 1].BeamNumber = beamcount
            ds.BeamSequence[beamcount - 1].TreatmentDeliveryType = 'TREATMENT'
            ds.BeamSequence[beamcount - 1].ReferencedPatientSetupNumber = beamcount 
            ds.BeamSequence[beamcount - 1].SourceAxisDistance = '1000'
            ds.BeamSequence[beamcount - 1].FinalCumulativeMetersetWeight = '1'
            ds.BeamSequence[beamcount - 1].PrimaryDosimeterUnit = 'MU'
            ds.BeamSequence[beamcount - 1].PrimaryFluenceModeSequence = Sequence()
            PrimaryFluenceMode1 = Dataset()
            ds.BeamSequence[beamcount - 1].PrimaryFluenceModeSequence.append(PrimaryFluenceMode1)
            ds.BeamSequence[beamcount - 1].PrimaryFluenceModeSequence[0].FluenceMode = 'STANDARD'
        if "      Name =" == line[:12] and beginbeam and noname:
            ds.BeamSequence[beamcount - 1].BeamName = re.findall(r'"([^"]*)"', line)[0]
            noname = False
        if "   PrescriptionPointName" in line:
            nameofrefpt = re.findall(r'"([^"]*)"', line)[0]
            for i, name in enumerate(point_names, 0):
                if nameofrefpt == name:
                    doserefpt = point_values[i]
            if doserefpt != []:
                print("Dose reference point: " + str([float(doserefpt[0])-xshift, float(doserefpt[1])-yshift, float(doserefpt[2])]))
                ds.FractionGroupSequence[0].ReferencedBeamSequence[beamcount - 1].BeamDoseSpecificationPoint = [float(doserefpt[0])-xshift, float(doserefpt[1])-yshift, float(doserefpt[2])] #Not sure if I need shifts here or not...?
            else:
                print("No dose reference point, setting to isocenter: " + str([float(isocenter[0]) - xshift, float(isocenter[1]) - yshift, float(isocenter[2])]))
                ds.FractionGroupSequence[0].ReferencedBeamSequence[beamcount - 1].BeamDoseSpecificationPoint = [float(isocenter[0]) - xshift, float(isocenter[1]) - yshift, float(isocenter[2])]
        if "      PrescriptionDose =" == line[:24]:
            prescdose = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "      Modality =" == line[:16] and beginbeam:
            if (re.findall(r'"([^"]*)"', line)[0] == 'Photons'):
                ds.BeamSequence[beamcount - 1].RadiationType = 'PHOTON'
            elif (re.findall(r'"([^"]*)"', line)[0] == 'Electrons'):
                ds.BeamSequence[beamcount - 1].RadiationType = 'ELECTRON'
            else:
                ds.BeamSequence[beamcount - 1].RadiationType = ""
        if "      SetBeamType" in line and beginbeam:
            if "STATIC" ==  re.findall(r'"([^"]*)"', line)[0].upper():
                ds.BeamSequence[beamcount - 1].BeamType = re.findall(r'"([^"]*)"', line)[0].upper()
            else:
                if "Step & Shoot" in re.findall(r'"([^"]*)"', line)[0] or ("step" in re.findall(r'"([^"]*)"', line)[0] and 'shoot' in re.findall(r'"([^"]*)"', line)[0]) or ("Step" in re.findall(r'"([^"]*)"', line)[0] and 'Shoot' in re.findall(r'"([^"]*)"', line)[0]):
                    #ds.BeamSequence[beamcount - 1].BeamType = "STATIC"
                    flag_stepnshoot = True
                #else:
                ds.BeamSequence[beamcount - 1].BeamType = "DYNAMIC"
        if "MonitorUnitInfo ={" in line and beginbeam:
            MUlineflag = True
            ctrlptmeterflag = False
        if "SourceToPrescriptionPointDistance" in line and MUlineflag == True:
            sad = float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])*10
            #ds.BeamSequence[beamcount - 1].SourceAxisDistance = sad
        if "PrescriptionDose =" in line and MUlineflag == True:
            prescripdose =  float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])
            normdose = float(re.findall(r"[-+]?\d*\.\d+|\d+", all_lines[linenum + 15])[0])
            OFc = float(re.findall(r"[-+]?\d*\.\d+|\d+", all_lines[linenum + 17])[0])
            # OFc value added by Achraf Touzani 2018
            if normdose == 0:
                beammu = 0 
                continue
            ds.FractionGroupSequence[0].ReferencedBeamSequence[beamcount - 1].BeamDose = float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])/100
            if beamenergies[beamcount-1] == '6': 
                beammu = float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])/(normdose*PDD6MV*OFc)
            elif beamenergies[beamcount-1] == '15': 
                beammu = float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])/(normdose*PDD15MV*OFc)
            elif beamenergies[beamcount-1] == '16':
                beammu = float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])/(normdose*PDD16MV*OFc)
            elif beamenergies[beamcount-1] == '10':
                beammu = float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])/(normdose*PDD10MV*OFc)
            else:
                print("\n \n Error, beam energy not 6, 10, 15 or 16 MV")
                return
            print("Beam MU: " + str(beammu))
            ds.FractionGroupSequence[0].ReferencedBeamSequence[beamcount - 1].BeamMeterset = beammu
            beamdoses.append(beammu)

            MUlineflag = False
            #Figure out what to do with BeamDose
        if "MachineNameAndVersion =" in line and nomachinename:
            machinename = re.findall(r'"([^"]*)"', line)[0]
            machinename = machinename.partition(":")[0]
            ds.BeamSequence[beamcount - 1].TreatmentMachineName = machinename
            nomachinename = False
        if "MachineEnergyName =" in line and beginbeam:
            beamenergies.append(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])
        if "NumberOfControlPoints" in line and beginbeam:
            numctrlpts = int(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])
            ds.BeamSequence[beamcount - 1].ControlPointSequence = Sequence()
            currentmeterset = 0.0
        if "ControlPointList ={" in line:
            ctrlptlist = True
            #print("ctrlptlist is True")
        if "Gantry =" in line and ctrlptlist: #Find out if this is gantry angle.
            gantryangle = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
            currentcontrolpoint = currentcontrolpoint + 1
        if "  Collimator =" in line and ctrlptlist:
            colangle = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "  Couch =" in line and ctrlptlist:
            psupportangle = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
        if "     WedgeName = " in line and ctrlptlist:
            print("wedge name found")
            if re.findall(r'"([^"]*)"', line)[0] == 'No Wedge' or re.findall(r'"([^"]*)"', line)[0] == "":
                wedgeflag = False
                print("Wedge is no name")
                numwedges = 0 
            elif "edw" in re.findall(r'"([^"]*)"', line)[0] or "EDW" in re.findall(r'"([^"]*)"', line)[0]:
                print("Wedge present")
                wedgetype = "DYNAMIC"
                wedgeflag = True
                numwedges = 1
                wedgeangle = re.findall(r"[-+]?\d*\.\d+|\d+", all_lines[linenum + 4])[0]
                wedgeinorout = ""
                wedgeinorout = re.findall(r'"([^"]*)"', all_lines[linenum+1])[0]
                if "WedgeBottomToTop" == wedgeinorout:
                    wedgename = re.findall(r'"([^"]*)"', line)[0].upper() +  wedgeangle + "IN"
                    wedgeorientation = '0' # temporary until I find out what to put here
                elif "WedgeTopToBottom" == wedgeinorout:
                    wedgename = re.findall(r'"([^"]*)"', line)[0].upper() +  wedgeangle + "OUT"
                    wedgeorientation = '180'
                print("Wedge name = ", wedgename)
            elif "UP" in re.findall(r'"([^"]*)"', line)[0]:
                print("Wedge present")
                wedgetype = "STANDARD"
                wedgeflag = True
                numwedges = 1
                wedgeangle = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
                wedgeinorout = ""
                wedgeinorout = re.findall(r'"([^"]*)"', all_lines[linenum+1])[0]
                if int(wedgeangle) == 15:
                    numberinname = '30'
                elif int(wedgeangle) == 45:
                    numberinname = '20'
                elif int(wedgeangle) == 30:
                    numberinname = '30'
                elif int(wedgeangle) == 60:
                    numberinname = '15'
                if "WedgeRightToLeft" == wedgeinorout:
                    wedgename = "W" +  str(int(wedgeangle)) + "R" + numberinname# + "U"
                    wedgeorientation = '90' # temporary until I find out what to put here
                elif "WedgeLeftToRight" == wedgeinorout:
                    wedgename = "W" +  str(int(wedgeangle)) + "L" + numberinname# + "U"
                    wedgeorientation = '270'
                elif "WedgeTopToBottom" == wedgeinorout:
                    wedgename = "W" +  str(int(wedgeangle)) + "OUT" + numberinname# + "U"
                    wedgeorientation = '180' # temporary until I find out what to put here
                elif "WedgeBottomToTop" == wedgeinorout:
                    wedgename = "W" +  str(int(wedgeangle)) + "IN" + numberinname# + "U"
                    wedgeorientation = '0' # temporary until I find out what to put here
                print("Wedge name = ", wedgename)
        if "LeftJawPosition" in line and x1 == "":
            x1 = str(-float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])*10)
            #print("X jaw 1:", x1, "\n")
        if "RightJawPosition" in line and x2 == "":
            x2 = float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])*10
            #print("X jaw 2:", x2, "\n")
        if "TopJawPosition" in line and y2 == "":
            y2 = float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])*10
            #print("y jaw 2:", y2, "\n")
        if "BottomJawPosition" in line and y1 == "":
            y1 = -float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])*10
            #print("y jaw 1:", y1, "\n")
        if "MLCLeafPositions ={" in line:
            mlcleafpos = True
            ctrlptmeterflag = True
        if mlcleafpos and "Points[] ={" in line:
            beginleafpoints = True
            del leafpositions1
            del leafpositions2
            leafpositions1 = []
            leafpositions2 = []
            continue
        if beginleafpoints:
            countpoints = countpoints + 1
            leafpointline = line.strip()
            leafpoints = leafpointline.split(',')
            #print("leafpoints: ", leafpoints)
            if leafpoints[0] == '};':
                beginleafpoints = False
                mlcleafpos = False
                leafpositions1 = list(reversed(leafpositions1))
                leafpositions2 = list(reversed(leafpositions2))
                totalleafpositions.append((leafpositions1+leafpositions2))
                continue
            leafpositions1.append(-float(leafpoints[0])*10)
            leafpositions2.append(float(leafpoints[1])*10)
        if ctrlptmeterflag and "  Weight =" in line:
            metersetweight.append(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])
        if "SSD = " in line and beginbeam:
            ssd = float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0])*10
        if "GantryIsCCW = " in line:
            if re.findall(r"[-+]?\d*\.\d+|\d+", line)[0] == '0':  #This may be a problem here!!!! Not sure how to Pinnacle does this, could be 1 if CW, must be somewhere that states if gantry is rotating or not
                gantryrotdir = 'NONE'
            elif  re.findall(r"[-+]?\d*\.\d+|\d+", line)[0] == '1':
                gantryrotdir = 'CC'
        if "GantryIsCW = " in line:
            if re.findall(r"[-+]?\d*\.\d+|\d+", line)[0] == '0':
                gantryrotdir = 'NONE'
            elif  re.findall(r"[-+]?\d*\.\d+|\d+", line)[0] == '1':
                gantryrotdir = 'CW'
        if flag_stepnshoot and "      DisplayMAXLeafMotion" in line:
            doserate = "400" 
            ds.BeamSequence[beamcount - 1].NumberOfControlPoints = numctrlpts*2
            ds.BeamSequence[beamcount - 1].SourceToSurfaceDistance = ssd
            if numwedges > 0:
                ds.BeamSequence[beamcount - 1].WedgeSequence = Sequence()
                Wedge1 = Dataset()
                ds.BeamSequence[beamcount - 1].WedgeSequence.append(Wedge1)
                ds.BeamSequence[beamcount-1].WedgeSequence[0].WedgeNumber = 1 #I am assuming only one wedge per beam (which makes sense because you can't change it during beam)
                ds.BeamSequence[beamcount-1].WedgeSequence[0].WedgeType = wedgetype #might need to change this
                ds.BeamSequence[beamcount-1].WedgeSequence[0].WedgeAngle = wedgeangle
                ds.BeamSequence[beamcount-1].WedgeSequence[0].WedgeID = wedgename
                ds.BeamSequence[beamcount-1].WedgeSequence[0].WedgeOrientation = wedgeorientation
                ds.BeamSequence[beamcount-1].WedgeSequence[0].WedgeFactor = ""
            
            #ds.BeamSequence[beamcount - 1].SourceAxisDistance = '1000'
            metercount = 1
            for j in range(0,numctrlpts*2):
                exec("ControlPoint%d = Dataset()"% (j+1))
                exec("ds.BeamSequence[beamcount - 1].ControlPointSequence.append(ControlPoint%d)"%(j+1))
                ds.BeamSequence[beamcount - 1].ControlPointSequence[j].ControlPointIndex = j
                ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence = Sequence()
                ds.BeamSequence[beamcount - 1].ControlPointSequence[j].ReferencedDoseReferenceSequence = Sequence()
                ReferencedDoseReference1 = Dataset()
                ds.BeamSequence[beamcount - 1].ControlPointSequence[j].ReferencedDoseReferenceSequence.append(ReferencedDoseReference1)
                if j%2 == 1: # odd number control point
                    curretnmeterset = currentmeterset + float(metersetweight[metercount])
                    metercount = metercount + 1

                ds.BeamSequence[beamcount - 1].ControlPointSequence[j].CumulativeMetersetWeight = currentmeterset
                ds.BeamSequence[beamcount - 1].ControlPointSequence[j].ReferencedDoseReferenceSequence[0].CumulativeDoseReferenceCoefficient = currentmeterset
                ds.BeamSequence[beamcount - 1].ControlPointSequence[j].ReferencedDoseReferenceSequence[0].ReferencedDoseReferenceNumber = '1'
                
                if j == 0: #first control point beam meterset always zero
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].NominalBeamEnergy = beamenergies[beamcount - 1]
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].DoseRateSet = doserate
                    
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].GantryRotationDirection = 'NONE'
                    #print("Gantry angle list length: ", len(gantryangles))
                    #print("current controlpoint: ", j)
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].GantryAngle = gantryangle
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDeviceAngle = colangle
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDeviceRotationDirection = 'NONE'
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].SourceToSurfaceDistance = ssd
                    BeamLimitingDevicePosition1 = Dataset() #This will be the x jaws
                    BeamLimitingDevicePosition2 = Dataset() #this will be the y jaws
                    if numwedges > 0:
                        WedgePosition1 = Dataset()
                        ds.BeamSequence[beamcount - 1].ControlPointSequence[j].WedgePositionSequence = Sequence()
                        ds.BeamSequence[beamcount - 1].ControlPointSequence[j].WedgePositionSequence.append(WedgePosition1)
                        ds.BeamSequence[beamcount - 1].ControlPointSequence[j].WedgePositionSequence[0].WedgePosition = "IN"
                        ds.BeamSequence[beamcount - 1].ControlPointSequence[j].WedgePositionSequence[0].ReferencedWedgeNumber = '1'
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence.append(BeamLimitingDevicePosition1)
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence.append(BeamLimitingDevicePosition2)
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence[0].RTBeamLimitingDeviceType = 'ASYMX'
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence[0].LeafJawPositions = [x1,x2]
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence[1].RTBeamLimitingDeviceType = 'ASYMY'
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence[1].LeafJawPositions = [y1, y2]
                    BeamLimitingDevicePosition3 = Dataset() #this will be the MLC 
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence.append(BeamLimitingDevicePosition3)
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence[2].RTBeamLimitingDeviceType = 'MLCX'
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence[2].LeafJawPositions = totalleafpositions[j]
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].SourceToSurfaceDistance = ssd
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDeviceRotationDirection = 'NONE'
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].PatientSupportAngle = psupportangle
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].PatientSupportRotationDirection = 'NONE'
                    print("Setting Isocenter postion: " + "[" + str(float(isocenter[0]) - xshift) +" , " +str(float(isocenter[1]) - yshift) + " , " + str(float(isocenter[2]))+ "]")
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].IsocenterPosition = [float(isocenter[0]) - xshift, float(isocenter[1]) - yshift, float(isocenter[2])]
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].GantryRotationDirection = gantryrotdir
                else:
                    BeamLimitingDevicePosition1 = Dataset() #This will be the mlcs for control points other than the first
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence.append(BeamLimitingDevicePosition1)
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence[0].RTBeamLimitingDeviceType = 'MLCX'
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence[0].LeafJawPositions = totalleafpositions[int(j/2)]
                ds.BeamSequence[beamcount - 1].NumberOfWedges = numwedges
                ds.BeamSequence[beamcount - 1].NumberOfCompensators = '0' # this is temporary value, will read in from file later
                ds.BeamSequence[beamcount - 1].NumberOfBoli = '0' # Also temporary
                ds.BeamSequence[beamcount - 1].NumberOfBlocks = '0' # Temp 
                ds.BeamSequence[beamcount - 1].BeamLimitingDeviceSequence = Sequence()
                BeamLimitingDevice1 = Dataset()
                BeamLimitingDevice2 = Dataset()
                BeamLimitingDevice3 = Dataset()
                ds.BeamSequence[beamcount - 1].BeamLimitingDeviceSequence.append(BeamLimitingDevice1)
                ds.BeamSequence[beamcount - 1].BeamLimitingDeviceSequence.append(BeamLimitingDevice2)
                ds.BeamSequence[beamcount - 1].BeamLimitingDeviceSequence.append(BeamLimitingDevice3)
                ds.BeamSequence[beamcount - 1].BeamLimitingDeviceSequence[0].RTBeamLimitingDeviceType = 'ASYMX'
                ds.BeamSequence[beamcount - 1].BeamLimitingDeviceSequence[1].RTBeamLimitingDeviceType = 'ASYMY'
                ds.BeamSequence[beamcount - 1].BeamLimitingDeviceSequence[2].RTBeamLimitingDeviceType = 'MLCX'
                ds.BeamSequence[beamcount - 1].BeamLimitingDeviceSequence[0].NumberOfLeafJawPairs = '1'
                ds.BeamSequence[beamcount - 1].BeamLimitingDeviceSequence[1].NumberOfLeafJawPairs = '1'
                ds.BeamSequence[beamcount - 1].BeamLimitingDeviceSequence[2].NumberOfLeafJawPairs = '60'
                bounds = ['-200','-190','-180','-170','-160','-150','-140','-130','-120','-110','-100','-95','-90','-85','-80','-75','-70','-65','-60','-55','-50','-45','-40','-35','-30','-25','-20','-15','-10','-5','0','5','10','15','20','25','30','35','40','45','50','55','60','65','70','75','80','85','90','95','100','110','120','130','140','150','160','170','180','190','200']
                ds.BeamSequence[beamcount - 1].BeamLimitingDeviceSequence[2].LeafPositionBoundaries = bounds
            ctrlptlist = False 
            wedgeflag = False
            numwedges = 0
            beginbeam = False
        if "      DisplayMAXLeafMotion" in line and not flag_stepnshoot:
            #doserate = re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]
            doserate = "400" 
            ds.BeamSequence[beamcount - 1].NumberOfControlPoints = numctrlpts + 1
            ds.BeamSequence[beamcount - 1].SourceToSurfaceDistance = ssd
            if numwedges > 0:
                ds.BeamSequence[beamcount - 1].WedgeSequence = Sequence()
                Wedge1 = Dataset()
                ds.BeamSequence[beamcount - 1].WedgeSequence.append(Wedge1)
                ds.BeamSequence[beamcount-1].WedgeSequence[0].WedgeNumber = 1 #I am assuming only one wedge per beam (which makes sense because you can't change it during beam)
                ds.BeamSequence[beamcount-1].WedgeSequence[0].WedgeType = wedgetype #might need to change this
                ds.BeamSequence[beamcount-1].WedgeSequence[0].WedgeAngle = wedgeangle
                ds.BeamSequence[beamcount-1].WedgeSequence[0].WedgeID = wedgename
                ds.BeamSequence[beamcount-1].WedgeSequence[0].WedgeOrientation = wedgeorientation
                ds.BeamSequence[beamcount-1].WedgeSequence[0].WedgeFactor = ""
            for j in range(0,numctrlpts+1):
                exec("ControlPoint%d = Dataset()"% (j+1))
                exec("ds.BeamSequence[beamcount - 1].ControlPointSequence.append(ControlPoint%d)"%(j+1))
                ds.BeamSequence[beamcount - 1].ControlPointSequence[j].ControlPointIndex = j
                ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence = Sequence()
                ds.BeamSequence[beamcount - 1].ControlPointSequence[j].ReferencedDoseReferenceSequence = Sequence()
                ReferencedDoseReference1 = Dataset()
                ds.BeamSequence[beamcount - 1].ControlPointSequence[j].ReferencedDoseReferenceSequence.append(ReferencedDoseReference1)
                ds.BeamSequence[beamcount - 1].ControlPointSequence[j].CumulativeMetersetWeight = metersetweight[j]
                if j == 0: #first control point beam meterset always zero
                    BeamLimitingDevicePosition1 = Dataset() #This will be the x jaws
                    BeamLimitingDevicePosition2 = Dataset() #this will be the y jaws
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].NominalBeamEnergy = beamenergies[beamcount - 1]
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].DoseRateSet = doserate
                    
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].GantryRotationDirection = 'NONE'
                    #print("Gantry angle list length: ", len(gantryangles))
                    #print("current controlpoint: ", j)
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].GantryAngle = gantryangle
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDeviceAngle = colangle
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].SourceToSurfaceDistance = ssd
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].ReferencedDoseReferenceSequence[0].CumulativeDoseReferenceCoefficient = '0'
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].ReferencedDoseReferenceSequence[0].ReferencedDoseReferenceNumber = '1'
                    if numwedges > 0:
                        WedgePosition1 = Dataset()
                        ds.BeamSequence[beamcount - 1].ControlPointSequence[j].WedgePositionSequence = Sequence()
                        ds.BeamSequence[beamcount - 1].ControlPointSequence[j].WedgePositionSequence.append(WedgePosition1)
                        ds.BeamSequence[beamcount - 1].ControlPointSequence[j].WedgePositionSequence[0].WedgePosition = "IN"
                        ds.BeamSequence[beamcount - 1].ControlPointSequence[j].WedgePositionSequence[0].ReferencedWedgeNumber = '1'
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence.append(BeamLimitingDevicePosition1)
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence.append(BeamLimitingDevicePosition2)
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence[0].RTBeamLimitingDeviceType = 'ASYMX'
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence[0].LeafJawPositions = [x1,x2]
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence[1].RTBeamLimitingDeviceType = 'ASYMY'
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence[1].LeafJawPositions = [y1, y2]
                    BeamLimitingDevicePosition3 = Dataset() #this will be the MLC 
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence.append(BeamLimitingDevicePosition3)
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence[2].RTBeamLimitingDeviceType = 'MLCX'
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence[2].LeafJawPositions = totalleafpositions[j]
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].SourceToSurfaceDistance = ssd
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDeviceRotationDirection = 'NONE'
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].PatientSupportAngle = psupportangle
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].PatientSupportRotationDirection = 'NONE'
                    print("No step-and-shoot Setting Isocenter postion: " + "[" + str(float(isocenter[0]) - xshift) +" , " +str(float(isocenter[1]) - yshift) + " , " + str(float(isocenter[2]))+ "]")
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].IsocenterPosition = [float(isocenter[0]) - xshift, float(isocenter[1]) - yshift, float(isocenter[2])]
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].GantryRotationDirection = gantryrotdir
                    ds.BeamSequence[beamcount - 1].NumberOfWedges = numwedges
                    ds.BeamSequence[beamcount - 1].NumberOfCompensators = '0' # this is temporary value, will read in from file later
                    ds.BeamSequence[beamcount - 1].NumberOfBoli = '0' # Also temporary
                    ds.BeamSequence[beamcount - 1].NumberOfBlocks = '0' # Temp 
                else:
                    BeamLimitingDevicePosition1 = Dataset() #This will be the mlcs for control points other than the first
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence.append(BeamLimitingDevicePosition1)
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence[0].RTBeamLimitingDeviceType = 'MLCX'
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].BeamLimitingDevicePositionSequence[0].LeafJawPositions = totalleafpositions[j-1]
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].ReferencedDoseReferenceSequence[0].CumulativeDoseReferenceCoefficient = '1'
                    ds.BeamSequence[beamcount - 1].ControlPointSequence[j].ReferencedDoseReferenceSequence[0].ReferencedDoseReferenceNumber = '1'
                
                ds.BeamSequence[beamcount - 1].BeamLimitingDeviceSequence = Sequence()
                BeamLimitingDevice1 = Dataset()
                BeamLimitingDevice2 = Dataset()
                BeamLimitingDevice3 = Dataset()
                ds.BeamSequence[beamcount - 1].BeamLimitingDeviceSequence.append(BeamLimitingDevice1)
                ds.BeamSequence[beamcount - 1].BeamLimitingDeviceSequence.append(BeamLimitingDevice2)
                ds.BeamSequence[beamcount - 1].BeamLimitingDeviceSequence.append(BeamLimitingDevice3)
                ds.BeamSequence[beamcount - 1].BeamLimitingDeviceSequence[0].RTBeamLimitingDeviceType = 'ASYMX'
                ds.BeamSequence[beamcount - 1].BeamLimitingDeviceSequence[1].RTBeamLimitingDeviceType = 'ASYMY'
                ds.BeamSequence[beamcount - 1].BeamLimitingDeviceSequence[2].RTBeamLimitingDeviceType = 'MLCX'
                ds.BeamSequence[beamcount - 1].BeamLimitingDeviceSequence[0].NumberOfLeafJawPairs = '1'
                ds.BeamSequence[beamcount - 1].BeamLimitingDeviceSequence[1].NumberOfLeafJawPairs = '1'
                ds.BeamSequence[beamcount - 1].BeamLimitingDeviceSequence[2].NumberOfLeafJawPairs = '60'
                bounds = ['-200','-190','-180','-170','-160','-150','-140','-130','-120','-110','-100','-95','-90','-85','-80','-75','-70','-65','-60','-55','-50','-45','-40','-35','-30','-25','-20','-15','-10','-5','0','5','10','15','20','25','30','35','40','45','50','55','60','65','70','75','80','85','90','95','100','110','120','130','140','150','160','170','180','190','200']
                ds.BeamSequence[beamcount - 1].BeamLimitingDeviceSequence[2].LeafPositionBoundaries = bounds
            ctrlptlist = False 
            wedgeflag = False
            numwedges = 0
            beginbeam = False
    ds.FractionGroupSequence[0].FractionGroupNumber = 1
    ds.FractionGroupSequence[0].NumberOfFractionsPlanned = numfracs
    ds.FractionGroupSequence[0].NumberOfBeams = beamcount
    ds.FractionGroupSequence[0].NumberOfBrachyApplicationSetups = '0'
    summed_pixel_values = []
    flag_nobinaryfile = False
    for currentbeam in range(0,beamcount):
        exec("PatientSetup%d = Dataset()"%(currentbeam+1))
        exec("ds.PatientSetupSequence.append(PatientSetup%d)"%(currentbeam+1))
        ds.PatientSetupSequence[currentbeam].PatientPosition = patient_position #get this from patient setup file
        ds.PatientSetupSequence[currentbeam].PatientSetupNumber = (currentbeam + 1)
        
        temp_pixelvalues, doseds = creatertdose(plannumber, planfolder, currentbeam + 1, beamdosefiles[currentbeam], beamdoses[currentbeam], numfracs)
        doseds.file_meta.MediaStorageSOPInstanceUID = doseinstuid + "." + str(plannumber)
        if flag_nobinaryfile:
            continue
        else:
            if currentbeam == 0:
                summed_pixel_values = temp_pixelvalues
            else:
                for i in range(0,len(summed_pixel_values)):
                    summed_pixel_values[i] = summed_pixel_values[i] + temp_pixelvalues[i]
    

    if flag_nobinaryfile == False:
        print("Max pixel value: " + str(max(summed_pixel_values)))
        print("Min pixel value: " + str(min(summed_pixel_values)))
        scale = max(summed_pixel_values) /65530
        doseds.DoseGridScaling = scale
        #doseds.TransferSyntaxUID=GTransferSyntaxUID
        #print(doseds.TransferSyntaxUID)

        print("Dose grid scaling: " + str(scale))
        
        ofile = open(Outputf + "/%s/samplebinaryslicevalues.txt"%(patientfolder),'w')
        pixel_binary_block = bytes()
        currline = 0
        pixelvaluelist = []
        for pp, element in enumerate(summed_pixel_values, 0):
            if(pp > dosexdim*doseydim*10 and pp < dosexdim*doseydim*11):
                currline = currline + 1
                ofile.write(str(element)+ " ")
                if currline%dosexdim == 0:
                    currline = 0
                    ofile.write("\n")
            if scale != 0:
                element = round(element/scale)
            else:
                element = 0
            pixelvaluelist.append(element)
            #pixel_binary_block += struct.pack("I", element)
        pixel_binary_block = struct.pack('%si' % len(pixelvaluelist), *pixelvaluelist)
        doseds.PixelData = pixel_binary_block
        ofile.close()
        dosefilename="RD."+doseds.file_meta.MediaStorageSOPInstanceUID+".dcm"
        print("\n Creating Dose file named : %s \n"%(dosefilename))
        doseds.save_as(Outputf+"%s/%s"%(patientfolder,dosefilename))
    #ds.FractionGroupSequence[0].ReferencedDoseReferenceSequence = Sequence()
    #ReferencedDoseReference2 = Dataset()
    #ds.FractionGroupSequence[0].ReferencedDoseReferenceSequence.append(ReferencedDoseReference2)
    #ds.FractionGroupSequence[0].ReferencedDoseReferenceSequence[0].TargetPrescriptionDose = int(prescdose)/int(numfracs)
    return ds
####################################################################################################################################################
####################################################################################################################################################


####################################################################################################################################################
#  Function: creatertdose
#  Purpose: create rt dose data structure and fill it
#  Requirements: Needs plan number and beam number
####################################################################################################################################################
def creatertdose(plannumber, planfolder, beamnum, binarynum, beamdosevalue, numfracs):
    global patientname
    global plansopinstuid
    global dob
    global pid
    global patient_sex
    global plansopinstuid
    global study_time
    global study_date
    global StudyInstanceUID
    global model
    global physician
    global doseseriesuid
    global doseinstuid
    global FrameUID
    global dosexdim
    global doseydim
    global dosezdim
    global doseoriginx
    global doseoriginy
    global doseoriginz
    global pixspacingy
    global pixspacingx
    global pixspacingz
    global posrefind
    global image_orientation
    global flag_nobinaryfile
    #Image Position (Patient) seems off, so going to calculate shift assuming dose origin in center and I want outer edge
    ydoseshift = float(pixspacingy)*float(doseydim)
    zdoseshift = float(pixspacingz)*float(dosezdim)
    #xdoseshift = float(pixspacingx)*float(dosexdim)/2
     # Populate required values for file meta information
    file_meta = Dataset()
    file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.481.2' # RT Dose Storage
    file_meta.TransferSyntaxUID = GTransferSyntaxUID
    file_meta.MediaStorageSOPInstanceUID = doseinstuid + "." + str(plannumber) + str(beamnum)
    file_meta.ImplementationClassUID = gImplementationClassUID #this value remains static since implementation for creating file is the same
    # Create the FileDataset instance (initially no data elements, but file_meta supplied)
    RDfilename="RD."+file_meta.MediaStorageSOPInstanceUID+".dcm"
    print("Dose file name : " + RDfilename)

    ds = FileDataset(RDfilename, {}, file_meta=file_meta, preamble=b'\x00'*128)
    
    ds.SpecificCharacterSet = 'ISO_IR 100'
    ds.InstanceCreationDate = time.strftime("%Y%m%d")
    ds.InstanceCreationTime = time.strftime("%H%M%S")
    ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.481.2' # RT Dose Storage
    ds.SOPInstanceUID = doseinstuid + "." + str(plannumber) + str(beamnum)
    ds.StudyDate = study_date
    ds.StudyTime = study_time
    ds.AccessionNumber = ''
    ds.Modality = 'RTDOSE'
    ds.Manufacturer = Manufacturer
    ds.OperatorsName = ""
    ds.ManufacturersModelName = model
    ds.SoftwareVersions = ['u\'9.0']
    ds.PhysiciansOfRecord = physician
    ds.PatientsName = patientname
    ds.PatientsBirthDate = dob
    ds.PatientID = pid
    ds.PatientsSex = patient_sex
    ds.SliceThickness = pixspacingz #Get this value from images???
    ds.StudyInstanceUID = StudyInstanceUID
    ds.SeriesInstanceUID = doseseriesuid + "." + str(plannumber) + str(beamnum)
    ds.StudyID = sid
    if(patient_position == 'HFS'):
        ds.ImagePositionPatient = [float(doseoriginx), float(doseoriginy) - ydoseshift , float(doseoriginz) - zdoseshift]
    elif(patient_position == 'HFP'):
        ds.ImagePositionPatient = [float(doseoriginx), float(doseoriginy) + ydoseshift , float(doseoriginz) - zdoseshift]
    elif(patient_position == 'FFS'):
        ds.ImagePositionPatient = [float(doseoriginx), float(doseoriginy) - ydoseshift , float(doseoriginz) + zdoseshift]
    elif(patient_position == 'FFP'):
        ds.ImagePositionPatient = [float(doseoriginx), float(doseoriginy) + ydoseshift , float(doseoriginz) + zdoseshift]
    ds.ImageOrientationPatient = image_orientation
    ds.FrameOfReferenceUID = FrameUID
    ds.PositionReferenceIndicator = posrefind #From image files?
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = 'MONOCHROME2'
    
    ds.NumberOfFrames = int(dosezdim) # is this Z dimension???
    ds.Rows = int(doseydim) #Using y for Rows because that's what's in the exported dicom file for test patient
    ds.Columns = int(dosexdim) #similar to above, x for columns
    ds.PixelSpacing = [pixspacingx, pixspacingy]
    ds.BitsAllocated = 32 #????
    ds.BitsStored = 32 #???
    ds.HighBit = 31 #???
    ds.PixelRepresentation = 0
    ds.DoseUnits =  'GY' #'RELATIVE'#'GY'
    ds.DoseType = 'PHYSICAL'
    ds.DoseSummationType = 'PLAN'
    ds.ReferencedRTPlanSequence = Sequence()
    ReferencedRTPlan1 = Dataset()
    ds.ReferencedRTPlanSequence.append(ReferencedRTPlan1)
    ds.ReferencedRTPlanSequence[0].ReferencedSOPClassUID = '1.2.840.10008.5.1.4.1.1.481.5'
    ds.ReferencedRTPlanSequence[0].ReferencedSOPInstanceUID = plansopinstuid + "." + str(plannumber)
    ds.ReferencedRTPlanSequence[0].ReferencedFractionGroupSequence = Sequence()
    ReferencedFractionGroup1 = Dataset()
    ds.ReferencedRTPlanSequence[0].ReferencedFractionGroupSequence.append(ReferencedFractionGroup1)
    ds.ReferencedRTPlanSequence[0].ReferencedFractionGroupSequence[0].ReferencedBeamSequence = Sequence()
    ReferencedBeam1 = Dataset()
    ds.ReferencedRTPlanSequence[0].ReferencedFractionGroupSequence[0].ReferencedBeamSequence.append(ReferencedBeam1)
    ds.ReferencedRTPlanSequence[0].ReferencedFractionGroupSequence[0].ReferencedBeamSequence[0].ReferencedBeamNumber = beamnum
    ds.ReferencedRTPlanSequence[0].ReferencedFractionGroupSequence[0].ReferencedFractionGroupNumber = '1'
    ds.TissueHeterogeneityCorrection = 'IMAGE'
    frameoffsetvect = []
    for p in range(0, int(dosezdim)):
        frameoffsetvect.append(int(p*int(pixspacingz)))
    ds.GridFrameOffsetVector = frameoffsetvect
    pixeldatallist = []
    randomcounter = 0
    print("Binary file: " + "plan.Trial.binary.%s"%binarynum)
    if os.path.isfile("%s%s/%s/plan.Trial.binary.%s"%(Inputf, patientfolder, planfolder, binarynum)):
        with open("%s%s/%s/plan.Trial.binary.%s"%(Inputf,patientfolder, planfolder, binarynum), "rb") as binary_file:
            data_element = binary_file.read(4)
            while data_element:
                randomcounter = randomcounter + 1
                value = struct.unpack(">f", data_element)[0]
                #if randomcounter % 5000 == 0:
                    #print("Prescibed dose: " + str(beamdosevalue))
                    #print("Pixel value: " + str(value*(beamdosevalue/100)*float(numfracs)))
                #value = value*(beamdosevalue)
                value = float(numfracs)*value*beamdosevalue/100
                pixeldatallist.append(value)
                data_element = binary_file.read(4)
    else:
        flag_nobinaryfile = True
    print("Length of Pixel Data list: " + str(len(pixeldatallist)))
    print("Z dim: " + str(dosezdim) + "       X dim: " + str(dosexdim)  + "       Y dim: " + str(doseydim))
    if len(pixeldatallist) == 0: #if the binary file is empty, treat as if it does not exist
        flag_nobinaryfile = True
    main_pix_array = []
    if flag_nobinaryfile == False:
        #ptag = ds.data_element("GridFrameOffsetVector").tag
        ds.FrameIncrementPointer = ds.data_element("GridFrameOffsetVector").tag
        
        for h in range(0, dosezdim):
            pixelsforframe = []
            for k in range(0, dosexdim*doseydim):
                #if(k > 0 and k%dosexdim == 0):
                #    topval = pixelsforframe[k - 1]
                #    for j in range(0, (dosexdim - 1)):
                #        pixelsforframe[k-j-1] = pixelsforframe[k-j-2]
                #    pixelsforframe[k-10] = topval
                pixelsforframe.append(float(pixeldatallist[h*doseydim*dosexdim + k]))
            #main_pix_array.append(pixelsforframe)
            main_pix_array = main_pix_array + list(reversed(pixelsforframe))

        main_pix_array = list(reversed(main_pix_array))

        temp_beamds = ds

        scale = max(main_pix_array) /65530
        temp_beamds.DoseGridScaling = scale
        #temp_beamds.TransferSyntaxUID=GTransferSyntaxUID
        pixel_binary_block = bytes()
        currline = 0
        pixelvaluelist = []
        for pp, element in enumerate(main_pix_array, 0):
        #if(pp > dosexdim*doseydim*10 and pp < dosexdim*doseydim*11):
            #currline = currline + 1
            #ofile.write(str(element)+ " ")
            #if currline%dosexdim == 0:
                #currline = 0
                #ofile.write("\n")
            if scale != 0:
                element = round(element/scale)
            else:
                element = 0
            pixelvaluelist.append(element)
        pixel_binary_block = struct.pack('%si' % len(pixelvaluelist), *pixelvaluelist)
        temp_beamds.PixelData = pixel_binary_block
        #print("\n Creating Dose file named : %s \n"%(RDfilename))
        #temp_beamds.save_as(Outputfolder+"%s/%s"%(patientfolder,RDfilename))
    return main_pix_array, ds

####################################################################################################################################################
####################################################################################################################################################
