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

import os
import sys
import SimpleITK as sitk
import numpy as np
from common import SaveImage, SaveDicomImage, LoadBValueImages, LoadDicomImage

class BValueModel:
    def __init__(self):
        self._imagesByBValue = None
        self._bValueImage = None
        self._outputPath = ""
        self._targetBValue = -1.0
        self._seriesNumber = 13701
        self._compress = False
        self._scale = 1.0

    def Name(self):
        raise RuntimeError("Not implemented.")

    def Good(self):
        if self.GetTargetBValue() < 0.0 or self._imagesByBValue is None or len(self._imagesByBValue) == 0:
            return False

        refImage = self._imagesByBValue[self._MinBValue()]

        return all((refImage.GetSize() == self._imagesByBValue[bValue].GetSize() for bValue in self._imagesByBValue))

    # Returns True if option is supported
    def SaveADC(self):
        return False

    # Returns True if option is supported
    def SavePerfusion(self):
        return False

    # Returns True if option is supported
    def SaveKurtosis(self):
        return False

    def SetCompress(self, compress):
        self._compress = compress

    def GetCompress(self):
        return self._compress

    def SetOutputPath(self, outputPath):
        self._outputPath = outputPath

    def GetOutputPath(self):
        return self._outputPath

    def GetOutputPathWithPrefix(self, prefix):
        if self.GetOutputPath().lower().endswith(".nii.gz"):
            outputPathNoExt = self.GetOutputPath()[:-7]
            ext = ".nii.gz"
        else:
            outputPathNoExt, ext = os.path.splitext(self.GetOutputPath())

        if len(ext) == 0:
            if os.altsep is not None:
                return self.GetOutputPath().rstrip(os.sep + os.altsep) + prefix
            else:
                return self.GetOutputPath().rstrip(os.sep) + prefix

        return outputPathNoExt + prefix + ext
        
    def GetADCOutputPath(self):
        return self.GetOutputPathWithPrefix("_ADC")

    def GetPerfusionOutputPath(self):
        return self.GetOutputPathWithPrefix("_Perfusion")

    def GetKurtosisOutputPath(self):
        return self.GetOutputPathWithPrefix("_Kurtosis")

    def SetImages(self, imagesByBValue):
        self._imagesByBValue = imagesByBValue
        return self._imagesByBValue is not None and len(self._imagesByBValue) > 0

    def GetImages(self):
        return self._imagesByBValue

    # Returns True if option is supported
    def SetADCImage(self, inputADCImage):
        return False

    def GetADCImage(self):
        return None

    def SetTargetBValue(self, bValue):
        self._targetBValue = bValue

    def GetTargetBValue(self):
        return self._targetBValue

    def SetBValueScale(self, scale):
        self._scale = scale

    def GetBValueScale(self):
        return self._scale

    def SetSeriesNumber(self, seriesNumber):
        self._seriesNumber = seriesNumber

    def GetSeriesNumber(self):
        return self._seriesNumber

    def GetADCSeriesNumber(self):
        return self.GetSeriesNumber() + 1

    def GetPerfusionSeriesNumber(self):
        return self.GetSeriesNumber() + 2

    def GetKurtosisSeriesNumber(self):
        return self.GetSeriesNumber() + 3

    def Run(self):
        raise RuntimeError("Not implemented.")

    def SaveImages(self):
        if self.GetOutputBValueImage() is None:
            return False

        print(f"Info: Saving b-value image to '{self.GetOutputPath()}' ...")

        if not self._SaveImage(self.GetOutputBValueImage(), self.GetOutputPath(), self.GetSeriesNumber(), f"{self.Name()}: Calculated b-{self.GetTargetBValue()}", self.GetTargetBValue()):
            print(f"Error: Failed to save b-value image.", file=sys.stderr)
            return False

        return True

    def GetOutputBValueImage(self):
        return self._bValueImage

    def _ComputeB0Image(self):
        # XXX: Call only this class' Good()
        if not BValueModel.Good(self):
            return False

        minBValue = self._MinBValue()

        if minBValue == 0.0:
            return True

        images = self.GetImages()

        model = MonoExponentialModel()

        model.SetImages(images)
        model.SetADCImage(self.GetADCImage())
        model.SetTargetBValue(0.0)

        print("Info: No b-0 image present. Precomputing b-0 image ...")

        if not model.Run():
            return False

        refImage = images[minBValue]

        b0Image = model.GetOutputBValueImage()
        b0Image.CopyInformation(refImage)

        for key in refImage.GetMetaDataKeys():
            b0Image.SetMetaData(key, refImage.GetMetaData(key))

        b0Image.SetMetaData("0018|9087", "0")
        images[0.0] = b0Image

        self.SetImages(images)

        print("Info: Done.")

        return self._MinBValue() == 0.0
        
    def _SaveImage(self, image, path, seriesNumber, seriesDescription, bValue = None):
        ext = os.path.splitext(path)[1]

        if len(ext) > 0:
            return SaveImage(image, path, self.GetCompress())

        try:
            tags = next((self.GetImages()[bValue] for bValue in self.GetImages() if self.GetImages()[bValue].HasMetaDataKey("0020|000e")))
        except:
            print(f"Error: Cannot save DICOM from non-DICOM images.", file=sys.stderr)
            return False

        if image.GetPixelID() == sitk.sitkFloat32 or image.GetPixelID() == sitk.sitkFloat64:
            npNewImage = np.clip(sitk.GetArrayViewFromImage(image)*1e3, -32768, 32767).astype(np.int16)
            newImage = sitk.GetImageFromArray(npNewImage)
            newImage.CopyInformation(image)
            image = newImage

        for key in tags.GetMetaDataKeys():
            image.SetMetaData(key, tags[key])

        if bValue is None:
            image.EraseMetaData("0018|9087")
        else:
            image.SetMetaData("0018|9087", str(bValue))

        image.EraseMetaData("0028|1050") # Window center
        image.EraseMetaData("0028|1051") # Window width
        image.EraseMetaData("0028|1055") # Window explanation

        # XXX: Prevent ITK from inverting a transform when saving... potentially causing loss of information in pixel intensity!
        image.SetMetaData("0028|1052", "0") # Rescale intercept
        image.SetMetaData("0028|1053", "1") # Rescale slope
        image.SetMetaData("0028|1054", "US") # US - "unspecified"

        image.SetMetaData("0028|0011", str(seriesNumber))
        image.SetMetaData("0008|103e", seriesDescription)

        # Derivative description
        image.SetMetaData("0008|2111", "PyComputeBValue")

        return SaveDicomImage(image, path, self.GetCompress())

    def _SaveADCImage(self, image, path, seriesNumber, seriesDescription):
        npNewImage = np.clip(sitk.GetArrayViewFromImage(image)*1e6, 0, 4095).astype(np.int16)
        newImage = sitk.GetImageFromArray(npNewImage)
        newImage.CopyInformation(image)

        for key in image.GetMetaDataKeys():
            newImage.SetMetaData(key, image.GetMetaData(key))

        return self._SaveImage(newImage, path, seriesNumber, seriesDescription)

    def _MinBValue(self):
        return min(self.GetImages().keys())

    def _MaxBValue(self):
        return max(self.GetImages().keys())

    def _GetBValueImage(self, bValue):
        if self.GetImages() is None:
            return None

        # Python's pitiful equivalent for std::find/std::find_if
        try:
            return next((self.GetImages()[b] for b in self.GetImages().keys() if b == bValue))
        except:
            return None

    def _GetBValuesSorted(self):
        if self.GetImages() is None:
            return []

        bValues = list(self.GetImages().keys())
        bValues.sort()

        return bValues

    def _GetLogIntensityImages(self):
        if self.GetImages() is None:
            return None

        minBValue = self._MinBValue()

        b0Image = self.GetImages()[minBValue]
        npB0Image = sitk.GetArrayViewFromImage(b0Image).astype(np.float32)

        npLogImages = np.zeros([len(self.GetImages()), npB0Image.size], dtype=np.float32)

        m = 0
        for bValue in self._GetBValuesSorted():
            if bValue == minBValue:
                continue

            bValueImage = self.GetImages()[bValue]
            npBValueImage = sitk.GetArrayFromImage(bValueImage).astype(np.float32)
            npMask = np.logical_and((npB0Image > 0), (npBValueImage < npB0Image))

            npBValueImage[npMask] = np.log(npBValueImage[npMask] / npB0Image[npMask])
            npBValueImage[np.logical_not(npMask)] = 0

            npLogImages[m,:] = npBValueImage.ravel()

            m += 1

        return npLogImages

class MonoExponentialModel(BValueModel):
    def __init__(self):
        super(MonoExponentialModel, self).__init__()

        self._saveADC = False
        self._adcImage = None
        self._inputADCImage = None

    def Name(self):
        return "Mono Exponential"

    def Good(self):
        if not super(MonoExponentialModel, self).Good():
            return False

        if self._inputADCImage is None:
            return len(self.GetImages()) > 1

        refImage = self.GetImages()[self._MinBValue()]

        return self._inputADCImage.GetSize() == refImage.GetSize()

    def SaveADC(self):
        self._saveADC = True
        return True

    def SetADCImage(self, inputADCImage):
        self._inputADCImage = inputADCImage
        return True

    def GetADCImage(self):
        return self._inputADCImage

    def Run(self):
        self._bValueImage = None
        self._adcImage = None

        if not self.Good():
            return False
  
        minBValue = self._MinBValue()
        targetBValue = self.GetTargetBValue()

        b0Image = self.GetImages()[minBValue]
        npB0Image = sitk.GetArrayViewFromImage(b0Image)

        self._bValueImage = self._GetBValueImage(targetBValue)

        solveADC = True
        solveBValue = True

        if self.GetADCImage() is not None:
            npADCImage = (sitk.GetArrayViewFromImage(self.GetADCImage()) / 1.0e6).astype(np.float32)
            self._adcImage = sitk.GetImageFromArray(npADCImage)

            self._adcImage.CopyInformation(self.GetADCImage())

            solveADC = False

        if self._bValueImage is not None:
            if self.GetBValueScale() != 1.0:
                npBValueImage = np.clip(self.GetBValueScale() * sitk.GetArrayViewFromImage(self._bValueImage), 0, 4095).astype(np.int16)
                self._bValueImage = sitk.GetImageFromArray(npBValueImage)
                self._bValueImage.CopyInformation(b0Image)

            solveBValue = False

        if not solveADC and not solveBValue:
            # Nothing to do
            return True
        elif not solveADC:
            npADCImage = sitk.GetArrayViewFromImage(self._adcImage)
            npBValueImage = np.clip(self.GetBValueScale() * npB0Image * np.exp(-(targetBValue - minBValue)*npADCImage), 0, 4095).astype(np.int16)
            self._bValueImage = sitk.GetImageFromArray(npBValueImage)

            self._bValueImage.CopyInformation(b0Image)
        elif not solveBValue:
            M = len(self.GetImages())-1
            N = 1
            K = npB0Image.size

            A = np.zeros([M,N], dtype=np.float32)
            B = -self._GetLogIntensityImages()

            # Possibly remove one row
            B = B[:M,:]

            m = 0
            for bValue in self._GetBValuesSorted():
                if bValue == minBValue:
                    continue

                A[m, 0] = bValue - minBValue
                m += 1

            X, res, rank, S = np.linalg.lstsq(A, B, rcond=None)

            npADCImage = np.reshape(X[0,:], npB0Image.shape)
            
            self._adcImage = sitk.GetImageFromArray(npADCImage)
            self._adcImage.CopyInformation(b0Image)
        else:
            M = len(self.GetImages())
            N = 2
            K = npB0Image.size

            A = np.zeros([M,N], dtype=np.float32)
            B = -self._GetLogIntensityImages()

            m = 0
            for bValue in self._GetBValuesSorted():
                if bValue == minBValue:
                    continue

                A[m, 0] = bValue - minBValue
                m += 1

            A[-1, 0] = targetBValue - minBValue
            A[-1, 1] = 1

            X, res, rank, S = np.linalg.lstsq(A, B, rcond=None)

            npADCImage = np.reshape(X[0,:], npB0Image.shape)
            npBValueImage = np.reshape(X[1,:], npB0Image.shape)

            npBValueImage = np.clip(self.GetBValueScale() * npB0Image * np.exp(npBValueImage), 0, 4095).astype(np.int16)

            self._adcImage = sitk.GetImageFromArray(npADCImage)
            self._bValueImage = sitk.GetImageFromArray(npBValueImage)

            self._adcImage.CopyInformation(b0Image)
            self._bValueImage.CopyInformation(b0Image)

        return True

    def SaveImages(self):
        if not super(MonoExponentialModel, self).SaveImages():
            return False

        if self._adcImage is not None and self._saveADC:
            print(f"Info: Saving ADC image to '{self.GetADCOutputPath()}' ...")

            if not self._SaveADCImage(self._adcImage, self.GetADCOutputPath(), self.GetADCSeriesNumber(), f"{self.Name()}: Calculated ADC"):
                print(f"Error: Failed to save ADC image.", file=sys.stderr)
                return False

        return True

