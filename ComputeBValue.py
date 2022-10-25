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
import argparse
from common import LoadBValueImages, ResolveBValueImages, LoadDicomImage, LoadImage
from models import MonoExponentialModel

modelTable = { 
  "mono": MonoExponentialModel
}

def main(modelType, outputPath, imagePaths, targetBValue, scale=1.0, seriesNumber=13701, saveADC=False, saveKurtosis=False, savePerfusion=False, compress=False, adcPath=None, initialBValue=0.0):
    if modelType not in modelTable:
        print(f"Error: Unknown model type '{modelType}'.", file=sys.stderr)
        exit(1)

    model = modelTable[modelType]()

    adcImage = None

    if adcPath is not None:
        if os.path.isdir(adcPath):
            adcImage = LoadDicomImage(adcPath)
        else:
            adcImage = LoadImage(adcPath)

        if adcImage is None:
            print(f"Error: Could not load ADC image '{adcPath}'.", file=sys.stderr)
            exit(1)

        print("Info: Loaded ADC image.")

        if not model.SetADCImage(adcImage):
            print(f"Warning: '{modelType}' model does not support using existing ADC image.", file=sys.stderr)

    images = []

    for imagePath in imagePaths:
        tmpImages = LoadBValueImages(imagePath)

        if tmpImages is None:
            print(f"Error: Could not load b-value image from '{imagePath}'.", file=sys.stderr)
            exit(1)

        images += tmpImages

    imagesByBValue = ResolveBValueImages(images, adcImage, initialBValue=initialBValue)

    if imagesByBValue is None:
        print(f"Error: Could not resolve b-values.", file=sys.stderr)
        exit(1)

    loadedBValues = list(imagesByBValue.keys())
    loadedBValues.sort()

    for bValue in loadedBValues:
        print(f"Info: Loaded b = {bValue}")

    model.SetTargetBValue(targetBValue)
    model.SetImages(imagesByBValue)
    model.SetOutputPath(outputPath)
    model.SetSeriesNumber(seriesNumber)
    model.SetCompress(compress)
    model.SetBValueScale(scale)

    if saveADC and not model.SaveADC():
        print(f"Warning: '{modelType}' model does not support saving ADC images.", file=sys.stderr)

    if savePerfusion and not model.SavePerfusion():
        print(f"Warning: '{modelType}' model does not support saving perfusion fraction images.", file=sys.stderr)

    if saveKurtosis and not model.SaveKurtosis():
        print(f"Warning: '{modelType}' model does not support saving kurtosis images.", file=sys.stderr)

    print(f"Info: Calculating b-{targetBValue}")

    if not model.Run():
        print("Error: Failed to compute b-value image.", file=sys.stderr)
        exit(1)

    if not model.SaveImages():
        print("Error: Failed to save output images.", file=sys.stderr)
        exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PyComputeBValue")
    parser.add_argument("-a", "--save-adc", dest="saveADC", action="store_true", default=False, help="Save calculated ADC. The output path will have _ADC appended (folder --> folder_ADC or file.ext --> file_ADC.ext).")
    parser.add_argument("-b", "--target-b-value", dest="targetBValue", required=True, type=float, help="Target b-value to calculate.")
    parser.add_argument("-c", "--compress", dest="compress", action="store_true", default=False, help="Compress output.")
    parser.add_argument("-k", "--save-kurtosis", dest="saveKurtosis", action="store_true", default=False, help="Save calculated kurtosis image. The output path will have _Kurtosis appended.")
    parser.add_argument("-n", "--series-number", dest="seriesNumber", type=int, default=13701, help="Series number for calculated b-value image.")
    parser.add_argument("-o", "--output-path", dest="outputPath", required=True, type=str, help="Output path which may be a folder for DICOM output or a medical image format file.")
    parser.add_argument("-p", "--save-perfusion", dest="savePerfusion", action="store_true", default=False, help="Save calculated perfusion fraction image. The output path will have _Perfusion appended.")
    parser.add_argument("-s", "--scale", dest="scale", type=float, default=1.0, help="Scale factor of target b-value image intensities.")
    parser.add_argument("-A", "--adc-path", dest="adcPath", required=False, type=str, default=None, help="Load an existing ADC image to use for computing a b-value image.")
    parser.add_argument("-I", "--iniital-b-value", dest="initialBValue", required=False, type=float, default=0.0, help="Initial expected b-value in a diffusion series of unknown b-values.")
    parser.add_argument("modelType", type=str, choices=list(modelTable.keys()), help="Diffusion model to use.")
    parser.add_argument("imagePaths", type=str, nargs='+', help="B-value diffusion series folders and image paths. Image paths may optionally be suffixed with ':bvalue' to indicate the diffusion b-value of the image. DICOM paths suffixed with ':-1' indicate that DICOM should be ignored when querying the b-value of the image.")

    args = parser.parse_args()

    main(**vars(args))


