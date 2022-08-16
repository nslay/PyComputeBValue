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
        for key in reader3D.GetMetaDataKeys(0): # Was 1
            image.SetMetaData(key, reader3D.GetMetaData(0, key)) # Was (1, key)

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

    manufacturer = manufacturer.strip().lower()

    if "siemens" in manufacturer:
        return _GetDiffusionBValueSiemens(img)
    elif "ge" in manufacturer:
        return _GetDiffusionBValueGE(img)
    elif "philips" in manufacturer:
        return _GetDiffusionBValuePhilips(img)

    return -1.0

def _Uninvert(img):
    pixelRep = _GetMetaData(img, "0028|0103")
    bitsStored = _GetMetaData(img, "0028|0101")

    if pixelRep is None or bitsStored is None:
        return img

    pixelRep = int(pixelRep)
    bitsStored = int(bitsStored)

    if bitsStored < 1:
        return img # Uhh?

    if pixelRep == 0: # Unsigned
        maxValue = (1 << bitsStored) - 1
    elif pixelRep == 1: # Signed
        maxValue = (1 << (bitsStored-1)) - 1
    else:
        return img # Uhh?

    halfValue = maxValue // 2

    npImg = sitk.GetArrayViewFromImage(img)

    # Not a b-value image or not "inverted" in the expected way
    if np.any(np.logical_or(npImg < 0, npImg > maxValue)):
        return img

    highCount = (npImg > halfValue).sum()
    #print(f"highCount = {highCount}, halfValue = {halfValue}, npImg.size = {npImg.size}")

    if 4*highCount > 3*npImg.size:
        print(f"Info: B-value image appears inverted. Restoring (high = {maxValue}) ...")
        npImg = maxValue - npImg

        newImg = sitk.GetImageFromArray(npImg)
        newImg.CopyInformation(img)

        for key in img.GetMetaDataKeys():
            newImg.SetMetaData(key, img.GetMetaData(key))

        return newImg

    return img

    
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

        if img is None:
            return None

        # Override bvalue
        img.SetMetaData("0018|9087", str(bValue))

        return [ _Uninvert(img) ]

    filesByBValue = ComputeBValueFileNames(path, seriesUID)        

    if filesByBValue is not None and len(filesByBValue) > 0:
        reader = sitk.ImageSeriesReader()
        reader.SetImageIO("GDCMImageIO")
        reader.SetLoadPrivateTags(True)
        reader.SetMetaDataDictionaryArrayUpdate(True)

        if dtype is not None:
            reader.SetOutputPixelType(dtype)

        images = []

        for bValue in filesByBValue:
            reader.SetFileNames(filesByBValue[bValue])

            try:
                img = reader.Execute()
            except:
                return None

            if not img.HasMetaDataKey("0020|000e"):
                for key in reader.GetMetaDataKeys(0):
                    img.SetMetaData(key, reader.GetMetaData(0, key))

            img.SetMetaData("0018|9087", str(bValue))

            images.append(_Uninvert(img))

        return images

    filesSortedByBValue = ComputeUnknownBValueFileNames(path, seriesUID)

    if filesSortedByBValue is not None and len(filesSortedByBValue) > 0:
        reader = sitk.ImageSeriesReader()
        reader.SetImageIO("GDCMImageIO")
        reader.SetLoadPrivateTags(True)
        reader.SetMetaDataDictionaryArrayUpdate(True)

        if dtype is not None:
            reader.SetOutputPixelType(dtype)

        images = []

        for imageFiles in filesSortedByBValue:
            reader.SetFileNames(imageFiles)

            try:
                img = reader.Execute()
            except:
                return None

            if not img.HasMetaDataKey("0020|000e"):
                for key in reader.GetMetaDataKeys(0):
                    img.SetMetaData(key, reader.GetMetaData(0, key))

            images.append(_Uninvert(img))

        return images

    return None

def ResolveBValueImages(images, adcImage, initialBValue = 0.0):
    if images is None or len(images) == 0:
        return None

    if all(image.HasMetaDataKey("0018|9087") for image in images):
        return { float(image.GetMetaData("0018|9087")): image for image in images }

    print(f"Info: Trying to infer unknown b-values (initial b = {initialBValue}) ...")

    if len(images) < 2:
        print("Error: Need at least two b-value images.", file=sys.stderr)
        return None

    if adcImage is None:
        print("Error: Need ADC image to infer unknown b-values.", file=sys.stderr)
        return None

    if any(image.GetSize() != adcImage.GetSize() for image in images):
        print("Error: Dimension mismatch between unknown b-value image and ADC.")
        return None

    lamIntensityKey = lambda image : np.percentile(sitk.GetArrayViewFromImage(image), 90)

    # Sort all images by intensity
    images.sort(key=lamIntensityKey, reverse=True)

    # Check that they're still sorted correctly
    prevBValue = -1.0
    for image in images:
        if image.HasMetaDataKey("0018|9087"):
            bValue = float(image.GetMetaData("0018|9087"))

            if bValue < prevBValue:
                print("Error: Known b-value image intensities are not consistent with b-value ordering.", file=sys.stderr)
                return None

            prevBValue = bValue

    # If the first b-value image has a b-value, use that instead of initialBValue
    if images[0].HasMetaDataKey("0018|9087"):
        initialBValue = float(images[0].GetMetaData("0018|9087"))
    else:
        images[0].SetMetaData("0018|9087", str(initialBValue))

    # Guess only first image was missing b-value?
    if all(image.HasMetaDataKey("0018|9087") for image in images):
        return { float(image.GetMetaData("0018|9087")): image for image in images }

    unknownIndices = [ i for i, image in enumerate(images) if not image.HasMetaDataKey("0018|9087") ]
    unknownIndices = { i: j for j, i in enumerate(unknownIndices) }

    numUnknowns = len(unknownIndices)
    numKnowns = len(images) - numUnknowns

    M = np.zeros([numUnknowns*numKnowns + numUnknowns*(numUnknowns-1)//2, numUnknowns])

    # Form the system
    r = 0
    for i in range(len(images)):
        for j in range(i+1, len(images)):
            if i not in unknownIndices and j not in unknownIndices:
                continue # Both known, no equation to solve

            if i in unknownIndices:
                M[r, unknownIndices[i]] = -1.0

            if j in unknownIndices:
                M[r, unknownIndices[j]] = 1.0

            r += 1

    logS = np.zeros(M.shape[:1])

    npADCImage = sitk.GetArrayViewFromImage(adcImage)/1.0e6

    # Form RHS
    r = 0
    for i in range(len(images)):
        bValueI = 0.0

        if images[i].HasMetaDataKey("0018|9087"):
            bValueI = float(images[i].GetMetaData("0018|9087"))

        npImageI = sitk.GetArrayViewFromImage(images[i])

        for j in range(i+1, len(images)):
            bValueJ = 0.0

            if images[j].HasMetaDataKey("0018|9087"):
                bValueJ = float(images[j].GetMetaData("0018|9087"))

            npImageJ = sitk.GetArrayViewFromImage(images[j])

            npMask = np.logical_and(npADCImage > 0, npImageJ > 0)
            npMask = np.logical_and(npMask, npImageI > 0)
            #npMask = np.logical_and(npMask, npImageI >= npImageJ)

            if not np.any(npMask):
                print("Error: Not enough valid voxels to solve for unknown b-value.", file=sys.stderr)
                return None

            # NOTE: bValueI and bValueJ are each 0 if they are unknown
            logMean = (npADCImage[npMask]*(npADCImage[npMask]*(bValueI - bValueJ) - np.log(npImageJ[npMask]/npImageI[npMask]))).mean()
            adcMean2 = (npADCImage[npMask]**2).mean()

            logS[r] = logMean / adcMean2
            r += 1

    #logS = np.inner(M.T, logS)
    #M = np.matmul(M.T, M)

    bValues, res, rank, S = np.linalg.lstsq(M, logS, rcond=None)

    imagesByBValue = dict()

    for i, image in enumerate(images):
        if i in unknownIndices:
            bValue = bValues[unknownIndices[i]]
        else:
            bValue = float(image.GetMetaData("0018|9087"))

        if bValue < initialBValue:
            print(f"Error: Unexpected solution for unknown b-value: b = {bValue}", file=sys.stderr)
            return None

        image.SetMetaData("0018|9087", str(bValue))

        if bValue in imagesByBValue:
            print(f"Error: Duplicate solved b-value {bValue}.", file=sys.stderr)
            return None

        imagesByBValue[bValue] = image

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

def ComputeUnknownBValueFileNames(path, seriesUID = None):
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

    if seriesUID is None or seriesUID == "":
        allSeriesUIDs = sitk.ImageSeriesReader.GetGDCMSeriesIDs(path)

        if len(allSeriesUIDs) == 0:
            return None

        seriesUID = allSeriesUIDs[0]

    fileNames = sitk.ImageSeriesReader.GetGDCMSeriesFileNames(path, seriesUID)

    if len(fileNames) == 0:
        return None

    filesAndImages = []

    for fileName in fileNames:
        try:
            reader.SetFileName(fileName)
            img = reader.Execute()
        except:
            return None

        filesAndImages.append((fileName, img))

    R = _GetOrientationMatrix(reader)

    if R is None:
        print(f"Error: Could not determine orientation matrix.", file=sys.stderr)
        return None

    # Sort slices by position
    # All the duplicate position slices will be contiguous
    lamPositionKey = lambda pair : np.inner(R.transpose(), _GetOrigin(pair[1]))[2]

    filesAndImages.sort(key=lamPositionKey)

    numBValues = 0
    i = 0
    inext = -1
    while i < len(filesAndImages):
        T = _GetOrigin(filesAndImages[i][1])

        # Identify contiguous range of slices with same position
        try:
            inext = next(iter(j for j, pair in enumerate(filesAndImages[i:], start=i) if np.linalg.norm(T - _GetOrigin(pair[1])) > 1e-10))
        except StopIteration:
            inext = len(filesAndImages)

        if numBValues == 0:
            numBValues = inext - i
        elif numBValues != inext - i:
            print(f"Error: Expected {numBValues} slices, but got {inext - i} (missing slice?).", file=sys.stderr)
            return None

        lamIntensityKey = lambda pair : np.percentile(sitk.GetArrayViewFromImage(pair[1]), 90)

        # Sort the slices by intensity

        filesAndImages[i:inext] = sorted(filesAndImages[i:inext], key=lamIntensityKey, reverse=True)

        i = inext

    filesSortedByBValue = []

    for i in range(numBValues):
        filesSortedByBValue.append([ pair[0] for pair in filesAndImages[i::numBValues] ])

    return filesSortedByBValue

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

