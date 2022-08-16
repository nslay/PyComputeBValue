# Introduction
PyComputeBValue is a tool that can calculate b-value images from one or more given b-value images (and optionally a given ADC image) using the mono-exponential model. It can read vendor-specific formatted diffusion b-values from DICOM, solve for unknown b-values, and can write b-value, Apparent Diffusion Coefficient (ADC), kurtosis, and perfusion fraction images to medical image formats like MetaIO, NIFTI, as well as DICOM. As diffusion MRI sequences are often interleaved with b-values, this tool supports extracting individual b-value image volumes from such sequences! When the tool cannot determine b-value from DICOM, it will attempt to solve for the unknown b-values (you must have ADC and know the initial b-value used, e.g. b=0).

PyComputeBValue is a Python port of ComputeBValue

https://github.com/nslay/ComputeBValue

NOTE: This tool is especially designed for prostate mpMRI imaging. Model assumptions follow those in the following paper:

Grant, Kinzya B., et al. "Comparison of calculated and acquired high b value diffusion-weighted imaging in prostate cancer." Abdominal imaging 40.3 (2015): 578-586.

# Installing
PyComputeBValue is a Python script and requires `SimpleITK` and `numpy`.

```
pip install SimpleITK numpy
```

# Usage
Clone the repository and run the script from the clone folder.

PyComputeBValue accepts one or more given DICOM folders, DICOM files or other medical images (e.g. NIFTI).

As a quickstart example

```
python ComputeBValue.py -b 1500 mono C:\Path\To\6-ep2ddifftraDYNDIST-03788
```

will compute a b-1500 image on ProstateX-0026's diffusion series using the mono exponential model. The output of this command might look like

```
Info: Loaded b = 50
Info: Loaded b = 400
Info: Loaded b = 800
Info: Calculating b-1500
Info: Saving b-value image to 'output.mha' ...
```

You may change the output path using the -o flag. Providing an output file path without a file extension will result in a DICOM folder. For example

```
python ComputeBValue.py -b 1500 -o B1500Image mono 6-ep2ddifftraDYNDIST-03788
```

will produce a DICOM folder with the following file hierarchy

```
B1500Image/
+-- 1.dcm
+-- 2.dcm
+-- 3.dcm
+-- 4.dcm
...
+-- 19.dcm
```

You may also manually provide the bvalue of any medical image by appending, ":bvalue", to the end of the path. For example:

```
python ComputeBValue.py -b 1500 -o B1500Image mono C:\Path\To\b50.nii.gz:50  C:\Path\To\b400.nii.gz:400 C:\Path\To\b800.nii.gz:800
```

This example reads B50, B400 and B800 images from NIFTI files.
  
This is useful for non-DICOM medical image formats that lack the extra meta-data as well as for DICOMs with missing information or undocumented vendor-specific tags for determining diffusion bvalue.

There are additional flags for saving the calculated ADC (-a), perfusion fraction (-p), and kurtosis image (-k). You may additionally compress the output images (-c) as well as change the output DICOM series number for the calculated b-value image (-n). By default, the series number is 13701 for the calculated b-value image. Other images share a similar series number offset by values of 1-3 (e.g. 13702-4).

**NOTE**: PyComputeBValue does not yet support any models that produce perfusion or kurtosis images.

As a special note, the perfusion fraction and kurtosis images are stored in non-DICOM medical images in floating point form. However, owing to limitations with ITK, perfusion and kurtosis image pixel values are scaled by 1000 and stored as integer values in DICOM. Both b-value and ADC images are unaffected and remain stored as integer values in both DICOM and non-DICOM formats.

Lastly, PyComputeBValue provides the below usage message when provided with the -h flag or no arguments. It's useful if you
forget.

```
usage: ComputeBValue.py [-h] [-a] -b TARGETBVALUE [-c] [-k] [-n SERIESNUMBER] -o OUTPUTPATH [-p] [-s SCALE] [-A ADCPATH] [-I INITIALBVALUE] {mono} imagePaths [imagePaths ...]

PyComputeBValue

positional arguments:
  {mono}                Diffusion model to use.
  imagePaths            B-value diffusion series folders and image paths. Image paths may optionally be suffixed with ':bvalue' to indicate the diffusion b-value of the image.

optional arguments:
  -h, --help            show this help message and exit
  -a, --save-adc        Save calculated ADC. The output path will have _ADC appended (folder --> folder_ADC or file.ext --> file_ADC.ext).
  -b TARGETBVALUE, --target-b-value TARGETBVALUE
                        Target b-value to calculate.
  -c, --compress        Compress output.
  -k, --save-kurtosis   Save calculated kurtosis image. The output path will have _Kurtosis appended.
  -n SERIESNUMBER, --series-number SERIESNUMBER
                        Series number for calculated b-value image.
  -o OUTPUTPATH, --output-path OUTPUTPATH
                        Output path which may be a folder for DICOM output or a medical image format file.
  -p, --save-perfusion  Save calculated perfusion fraction image. The output path will have _Perfusion appended.
  -s SCALE, --scale SCALE
                        Scale factor of target b-value image intensities.
  -A ADCPATH, --adc-path ADCPATH
                        Load an existing ADC image to use for computing a b-value image.
  -I INITIALBVALUE, --initial-b-value INITIALBVALUE
                        Initial expected b-value in a diffusion series of unknown b-values.
```

# Models
This section briefly describes how the models are implemented. All models solve a least-squares minimization problem of the form:

```
\ell(\theta) = \sum_{b \in B} (f(\theta) - \log(S_b/S_0))^2
```

where B is the set of b-values, `S_b` is a pixel intensity from a b-value
image, and `f(\theta)` is the model predicting exponential terms. In 
addition to `\theta`, one `S_b` may be an unknown value to optimize for.

The models are discussed in more detail in the following paper:

Grant, Kinzya B., et al. "Comparison of calculated and acquired high b value diffusion-weighted imaging in prostate cancer." Abdominal imaging 40.3 (2015): 578-586.

## Mono Exponential
The mono exponential model relates `S_b` and `D` through the expression:

```
S_b = S_0 e^{-bD}
```
where `D` is the ADC value.

This gives the parameters `\theta = \{ D \}` and the model function

```
f(\theta) = -bD
```

**NOTE**: This model can cope with the absence of `S_0` using a mathematical
trick. In this model, the ratio of two b-value image intensities is

```
(S_a / S_0) / (S_b / S_0) = e^{-(a-b)D}
```
With some rearrangement we have:

```
S_a = S_b e^{-(a-b)D}
```
Hence, by simply shifting by the minimum b-value, this model may be
used to calculate any b-value image in the absence of a proper b-0
image (it may even be used to calculate the b-0 image!).

IVIM and Kurtosis models may be ported in the future.
