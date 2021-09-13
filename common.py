# 
# Nathan Lay
# AI Resource at National Cancer Institute
# National Institutes of Health
# September 2021
# 
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR(S) ``AS IS'' AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
# NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 

import sys
import os
import re
import SimpleITK as sitk
import numpy as np
import uuid
import time
import base64

def SaveImage(image, path, compress=True):
    writer = sitk.ImageFileWriter()
    writer.SetFileName(path)
    writer.SetUseCompression(compress)

    try:
        writer.Execute(image)
    except:
        return False

    return True

def LoadImage(path, dtype = None):
    reader = sitk.ImageFileReader()
    reader.SetFileName(path)

    if dtype is not None:
        reader.SetOutputPixelType(dtype)

    try:
        return reader.Execute()
    except:
        return None

def LoadDicomImage(path, seriesUID = None, dim = None, dtype = None):
    if not os.path.exists(path):
        return None

    reader2D = sitk.ImageFileReader()
    reader2D.SetImageIO("GDCMImageIO")
    reader2D.SetLoadPrivateTags(True)

    if dtype is not None:
        reader2D.SetOutputPixelType(dtype)

    if dim is None: # Guess the dimension by the path
        dim = 2 if os.path.isfile(path) else 3

    if dim == 2:
        reader2D.SetFileName(path)

        try:
            return reader2D.Execute()
        except:
            return None

    if os.path.isfile(path):
        reader2D.SetFileName(path)

        try:
            reader2D.ReadImageInformation()
            seriesUID = reader2D.GetMetaData("0020|000e").strip()
        except:
            return None

        path = os.path.dirname(path)
        
    fileNames = []

    if seriesUID is None or seriesUID == "":
        allSeriesUIDs = sitk.ImageSeriesReader.GetGDCMSeriesIDs(path)

        if len(allSeriesUIDs) == 0:
            return None

        for tmpUID in allSeriesUIDs:
            tmpFileNames = sitk.ImageSeriesReader.GetGDCMSeriesFileNames(path, tmpUID)        
            
            if len(tmpFileNames) > len(fileNames):
                seriesUID = tmpUID
                fileNames = tmpFileNames # Take largest series
    else:
        fileNames = sitk.ImageSeriesReader.GetGDCMSeriesFileNames(path, seriesUID)

    if len(fileNames) == 0: # Huh?
        return None

    reader3D = sitk.ImageSeriesReader()
    reader3D.SetImageIO("GDCMImageIO")
    reader3D.SetFileNames(fileNames)
    reader3D.SetLoadPrivateTags(True)
    reader3D.SetMetaDataDictionaryArrayUpdate(True)

    #reader3D.SetOutputPixelType(sitk.sitkUInt16)

    if dtype is not None:
        reader3D.SetOutputPixelType(dtype)

    try:
        image = reader3D.Execute()
    except:
        return None

    # Check if meta data is available!
    # Copy it if it is not!
    if not image.HasMetaDataKey("0020|000e"):
        for key in reader3D.GetMetaDataKeys(1):
            image.SetMetaData(key, reader3D.GetMetaData(1, key))

    return image

def SaveDicomImage(image, path, compress=True):
    # Implement pydicom's behavior
    def GenerateUID(prefix="1.2.826.0.1.3680043.8.498."):
        if not prefix:
            prefix = "2.25."
    
        return str(prefix) + str(uuid.uuid4().int)

    if image.GetDimension() != 2 and image.GetDimension() != 3:
        raise RuntimeError("Only 2D or 3D images are supported.")

    if not image.HasMetaDataKey("0020|000e"):
        print("Error: Reference meta data does not appear to be DICOM?", file=sys.stderr)
        return False

    writer = sitk.ImageFileWriter()
    writer.SetImageIO("GDCMImageIO")
    writer.SetKeepOriginalImageUID(True)
    writer.SetUseCompression(compress)

    newSeriesUID = GenerateUID()

    if image.GetDimension() == 2:
        writer.SetFileName(path)

        imageSlice = sitk.Image([image.GetSize()[0], image.GetSize()[1], 1], image.GetPixelID(), image.GetNumberOfComponentsPerPixel())
        imageSlice.SetSpacing(image.GetSpacing())

        imageSlice[:,:,0] = image[:]

        # Copy meta data
        for key in image.GetMetaDataKeys():
            imageSlice.SetMetaData(key, image.GetMetaData(key))

        newSopInstanceUID = GenerateUID()

        imageSlice.SetMetaData("0020|000e", newSeriesUID)
        imageSlice.SetMetaData("0008|0018", newSopInstanceUID)
        imageSlice.SetMetaData("0008|0003", newSopInstanceUID)

        try:
            writer.Execute(image)
        except:
            return False

        return True

    if not os.path.exists(path):
        os.makedirs(path)

    for z in range(image.GetDepth()):
        newSopInstanceUID = GenerateUID()

        imageSlice = sitk.Image([image.GetSize()[0], image.GetSize()[1], 1], image.GetPixelID(), image.GetNumberOfComponentsPerPixel())

        imageSlice[:] = image[:,:,z]
        imageSlice.SetSpacing(image.GetSpacing())

        # Copy meta data
        for key in image.GetMetaDataKeys():
            imageSlice.SetMetaData(key, image.GetMetaData(key))

        # Then write new meta data ...
        imageSlice.SetMetaData("0020|000e", newSeriesUID)
        imageSlice.SetMetaData("0008|0018", newSopInstanceUID)
        imageSlice.SetMetaData("0008|0003", newSopInstanceUID)

        # Instance creation date and time
        imageSlice.SetMetaData("0008|0012", time.strftime("%Y%m%d"))
        imageSlice.SetMetaData("0008|0013", time.strftime("%H%M%S"))

        # Image number
        imageSlice.SetMetaData("0020|0013", str(z+1))

        position = image.TransformIndexToPhysicalPoint((0,0,z))

        # Image position patient
        imageSlice.SetMetaData("0020|0032", f"{position[0]}\\{position[1]}\\{position[2]}")

        # Slice location
        imageSlice.SetMetaData("0020|1041", str(position[2]))

        # Spacing
        imageSlice.SetMetaData("0018|0050", str(image.GetSpacing()[2]))
        imageSlice.SetMetaData("0018|0088", str(image.GetSpacing()[2]))

        imageSlice.EraseMetaData("0028|0106")
        imageSlice.EraseMetaData("0028|0107")

        slicePath = os.path.join(path, f"{z+1}.dcm")
        writer.SetFileName(slicePath)

        try:
            writer.Execute(imageSlice)
        except:
            print(f"Error: Failed to write slice '{slicePath}'.", file=sys.stderr)
            return False

    return True

def GetDiffusionBValue(img):
    tmp = _GetMetaData(img, "0018|9087")

    if tmp is not None:
        return float(tmp)

    # Check for ProstateX
    patientName = _GetMetaData(img, "0010|0010")
    patientId = _GetMetaData(img, "0010|0020")

    if patientName is None or patientId is None:
        return -1.0

    patientName = patientName.strip().lower()
    patientId = patientId.strip().lower()

    if "prostatex" in patientName or "prostatex" in patientId:
        return _GetDiffusionBValueProstateX(img)

    manufacturer = _GetMetaData(img, "0008|0070")

    if manufacturer is None:
        return -1.0

    manufacturer = manufacturer.strip.lower()

    if "siemens" in manufacturer:
        return _GetDiffusionBValueSiemens(img)
    elif "ge" in manufacturer:
        return _GetDiffusionBValueGE(img)
    elif "philips" in manufacturer:
        return _GetDiffusionBValuePhilips(img)

    return -1.0

def LoadBValueImages(path, seriesUID = None, dtype = None):
    # Check for hints
    tmp = re.search(":[0-9]+$", path)

    if tmp is not None:
        bValue = float(tmp.group(0)[1:])
        path = path[:tmp.start()]

        img = None
        if os.path.isdir(path):
            img = LoadDicomImage(path, seriesUID, dtype=dtype)
        else:
            img = LoadImage(path, dtype=dtype)

        # Override bvalue
        img.SetMetaData("0018|9087", str(bValue))

        return { bValue: img }

    filesByBValue = ComputeBValueFileNames(path, seriesUID)        

    if filesByBValue is None or len(filesByBValue) == 0:
        return None

    reader = sitk.ImageSeriesReader()
    reader.SetImageIO("GDCMImageIO")
    reader.SetLoadPrivateTags(True)
    reader.SetMetaDataDictionaryArrayUpdate(True)

    if dtype is not None:
        reader.SetOutputPixelType(dtype)

    imagesByBValue = dict()

    for bValue in filesByBValue:
        reader.SetFileNames(filesByBValue[bValue])

        try:
            img = reader.Execute()
        except:
            return None

        if not img.HasMetaDataKey("0020|000e"):
            for key in reader.GetMetaDataKeys(1):
                img.SetMetaData(key, reader.GetMetaData(1, key))

        imagesByBValue[bValue] = img

    return imagesByBValue

def ComputeBValueFileNames(path, seriesUID = None):
    reader = sitk.ImageFileReader()
    reader.SetImageIO("GDCMImageIO")
    reader.SetLoadPrivateTags(True)

    if not os.path.isdir(path):
        try:
            reader.SetFileName(path)
            reader.ReadImageInformation()

            seriesUID = reader.GetMetaData("0020|000e")
            path = os.path.dirname(path)
        except:
            return None

        if GetDiffusionBValue(reader) < 0.0:
            return None

    if seriesUID is None or seriesUID == "":
        allSeriesUIDs = sitk.ImageSeriesReader.GetGDCMSeriesIDs(path)

        if len(allSeriesUIDs) == 0:
            return None

        seriesUID = allSeriesUIDs[0]

    fileNames = sitk.ImageSeriesReader.GetGDCMSeriesFileNames(path, seriesUID)

    if len(fileNames) == 0:
        return None

    filesAndDictionariesByBValue = dict()

    for fileName in fileNames:
        try:
            reader.SetFileName(fileName)
            reader.ReadImageInformation()
        except:
            return None

        bValue = GetDiffusionBValue(reader)

        if bValue < 0.0:
            print(f"Error: Could not extract diffusion b-value for '{fileName}'.", file=sys.stderr)
            return None

        tags = { key: reader.GetMetaData(key) for key in reader.GetMetaDataKeys() }

        if bValue not in filesAndDictionariesByBValue:
            filesAndDictionariesByBValue[bValue] = [(fileName, tags)]
        else:
            filesAndDictionariesByBValue[bValue].append((fileName, tags))

    R = _GetOrientationMatrix(reader)

    if R is None:
        print(f"Error: Could not determine orientation matrix.", file=sys.stderr)
        return None

    # Sort slices by position
    lamPositionKey = lambda pair : np.inner(R.transpose(), _GetOrigin(pair[1]))[2]

    for bValue in filesAndDictionariesByBValue:
        filesAndDictionariesByBValue[bValue].sort(key=lamPositionKey)

    filesByBValue = dict()

    for bValue in filesAndDictionariesByBValue:
        filesByBValue[bValue] = [ pair[0] for pair in filesAndDictionariesByBValue[bValue] ]

    return filesByBValue

###### Private stuff below ######

# Find out more here:
# https://www.na-mic.org/wiki/NAMIC_Wiki:DTI:DICOM_for_DWI_and_DTI

def _GetMetaData(img, key):
    if isinstance(img, dict):
        return img[key] if key in img else None

    return img.GetMetaData(key) if img.HasMetaDataKey(key) else None

def _GetOrientationMatrix(img):
    tmp = _GetMetaData(img, "0020|0037")

    if tmp is None:
        return None

    R = np.zeros(9)
    R[0:6] = np.array([ float(x) for x in tmp.split("\\") ])
    R[6:9] = np.cross(R[0:3], R[3:6])

    R = np.reshape(R, [3,3])

    return R.transpose()


def _GetOrigin(img):
    tmp = _GetMetaData(img, "0020|0032")

    if tmp is None:
        return None

    return np.array([ float(x) for x in tmp.split("\\") ])

def _ExtractCSAHeader(img):
    return None

def _GetDiffusionBValueProstateX(img):
    sequenceName = _GetMetaData(img, "0018|0024")

    if sequenceName is None:
        return -1.0

    tmp = re.search("b[0-9]+t", sequenceName)

    if tmp is None:
        return -1.0

    return float(tmp.group(0)[1:-1])

def _GetDiffusionBValueSiemens(img):
    model = _GetMetaData(img, "0008|1090")

    if model is None:
        return -1.0

    model = model.strip().lower()
    
    if any((possibleModel in model for possibleModel in ["skyra", "verio"])):
        bvalue = _GetDiffusionBValueProstateX(img)

        if bvalue >= 0.0:
            return bvalue


    """
    csaDict = _ExtractCSAHeader(img)

    if not csaDict or "B_value" not in csaDict:
        return -1.0

    return float(csaDict["B_value"])
    """

    return -1.0

def _GetDiffusionBValueGE(img):
    tmp = _GetMetaData("0043|1039")

    if tmp is None:
        return -1.0

    tmp = tmp.split("\\")

    if len(tmp) <= 1:
        return -1.0 # What?

    try:
        return float(tmp[0]) % 100000
    except:
        return -1.0

    return -1.0 # Not reached

def _GetDiffusionBValuePhilips(img):
    tmp = _GetMetaData(img, "2001|1003")

    if tmp is None:
        return -1.0

    return float(tmp)

