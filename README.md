# MitoSLIT Pipeline

Welcome to **MitoSLIT** (*Mitochondria Segmentation with Locally Iterative Thresholding*)!

This pipeline is designed for researchers and students to segment mitochondrial networks from fluorescence microscopy images (such as confocal, TPEF, or 2p-FLIM TIFF images). The run time depends on the image size and resolution. The processing time will be longer if the image size is larger (more pixels) or the mitochondrial structures are more complicated. A typical running time for an image with 512*512 (pixels*pixels) is 40-80 s.

---

## Table of Contents
1. [Prerequisites & Installation](#1-prerequisites--installation)
2. [Quick Start (3 Steps)](#2-quick-start-3-steps)
3. [Image Organization & Pairing](#3-image-organization--pairing)
4. [User Configuration Guide](#4-user-configuration-guide)
5. [Understanding Output Files (4 Binary Masks)](#5-understanding-output-files-4-binary-masks)
6. [Speed vs. Precision Modes](#6-speed-vs-precision-modes)
7. [Troubleshooting & FAQs](#7-troubleshooting--faqs)

---

## 1. Prerequisites & Installation

MitoSLIT requires Python 3.10 to 3.12 and standard scientific image processing libraries.

### Option A: Using `pip` (Standard)
Install all required dependencies using [`requirements.txt`](./requirements.txt):
```bash
pip install -r requirements.txt
```

### Option B: Using Conda / Mamba (Recommended for Beginners)
Create and activate a dedicated environment using [`environment.yml`](./environment.yml):
```bash
conda env create -f environment.yml
conda activate mitoSLIT
```

---

## 2. Quick Start (3 Steps)

### Step 1: Put Your Images in a Folder
Place your raw `.tif` or `.tiff` images into an input folder (for example, `./MitoSLIT/`).

### Step 2: Configure the Script
Open [`MitoSLIT_pipeline.py`](./MitoSLIT_pipeline.py) in your code editor (e.g., VS Code, Spyder, PyCharm).  
Scroll to **Section 1: USER CONFIGURATION** (lines 85–148) and customize settings if needed:
- `input_dir`: Path to your raw images folder (e.g. `r".\MitoSLIT"`).
- `output_dir`: Path where output masks will be saved (e.g. `r".\MitoSLIT\filtered"`).
- `fluo_keyword`: Keyword that identifies your fluorescence images (e.g. `"FluoIntensity"` or `"FLIM"` or `"Ch2"` or `"GFP"`).
- `roi_keyword`: Keyword that identifies your tissue ROI masks (e.g. `"roi"`).

### Step 3: Run the Script
Run the pipeline from your terminal or command prompt:
```bash
python MitoSLIT_pipeline.py
```
You will see progress for each image in the terminal, and all segmented binary masks will be saved automatically to your `output_dir`.

---

## 3. Image Organization & Pairing

The pipeline automatically scans your `input_dir` and pairs each **fluorescence intensity image** with its corresponding **tissue ROI mask** based on your keywords.

### How Image Pairing Works
- **Fluorescence images**: Must contain `fluo_keyword` (case-insensitive) and must *not* contain `roi_keyword` or any `exclude_keywords`.
- **Tissue ROI images**: Must contain `roi_keyword` (case-insensitive).
- **Matching rule**: Matches filename stems. It also automatically handles slice numbers (e.g. `_01`, `_02` for z-stacks) and channel identifiers (e.g. `_Ch0`, `_Ch1`, `_Ch2`).

### Example File Naming:
| Fluorescence Image | Matching Tissue ROI Image | Result |
| :--- | :--- | :--- |
| `SAMPLE_WT_03_FluoIntensity.tif` | `SAMPLE_WT_03_roi.tif` | Paired automatically (all 4 masks saved) |
| `worm01_01_FLIM.tif` | `worm01_roi.tif` | Paired (handles z-slice `_01`) |
| `animal2_Ch2.tif` | `animal2_roi.tif` | Paired automatically |
| `sample_GFP.tif` | *(No ROI file)* | Processed normally (saves 2 unmasked images) |

> [!NOTE]
> In tissue ROI images, pixel value `0` represents the valid tissue area, while values `> 0` represent background/exterior regions that get masked out. If no matching ROI file is found for an image, the pipeline gracefully skips the ROI masking step and saves the 2 unmasked binary images.

---

## 4. User Configuration Guide

All user-configurable parameters are centralized in `MitoSLITConfig` near the top of [`MitoSLIT_pipeline.py`](./MitoSLIT_pipeline.py):

### File & Folder Paths
- `input_dir` *(str)*: Folder containing raw TIFF images. Default: `r".\MitoSLIT"`.
- `output_dir` *(str)*: Folder where segmented masks are saved. Default: `r".\MitoSLIT\filtered"`.

### Keyword Filters
- `fluo_keyword` *(str or tuple)*: Keyword(s) required in fluorescence filenames (case-insensitive).  
  *Examples*: `"FluoIntensity"`, `"FLIM"`, `"Ch2"`, `"GFP"`. Leave as `""` or `()` to match all TIFFs.
- `roi_keyword` *(str)*: Keyword identifying tissue ROI files (case-insensitive).  
  *Examples*: `"roi"`, `"_roi"`, `"TissueROI"`. Set to `""` if you do not have ROI files.
- `exclude_keywords` *(str or tuple)*: Keywords that ignore unwanted images (case-insensitive).  
  *Examples*: `("Ch1", "mask")` (ignores transmission channel `Ch1` or previously generated `mask` files).

### Denoising Method
- `denoise_method` *(int)*:
  - `0`: No denoising (uses raw image directly).
  - `1` **(Default & Recommended)**: Rolling-ball background subtraction + Gaussian smoothing. Ideal for images with uneven illumination or autofluorescence. Preserves camera integer depth to match `v4` bit-for-bit.
  - `2`: Grayscale morphological opening. Useful for removing high-frequency speckle noise.

### Core Segmentation Parameters
- `tile_size` *(int, pixels)*: Size of local analysis tiles. Should be larger than the largest single mitochondrial structure.  
  *Recommended*: `20` for somatic tissues (*C. elegans* intestine, muscle, aging); `10` for narrow structures like neurons.
- `smallest_mito_size` *(int, pixels)*: Minimum feature dimension. Objects smaller than `smallest_mito_size^2` pixels are removed as noise. Default: `2`.
- `local_window_scales` *(tuple of ints)*: Scale multipliers of `tile_size` for multiscale analysis. Default: `(1, 2)` (i.e. windows of 20 and 40 px).
- `skipping_pixels_list` *(tuple of ints)*: Sliding-window pixel skipping stride per scale.
  - `(0, 1)` **(Default)**: Standard mode matching `MitoSLIT_v4` 100% bit-for-bit (step sizes 1 and 2 px across windows).
  - `(1, 3)`: Fast mode (step sizes 2 and 4 px, >99.3% consistent, ~9s per image).
  - `(3, 7)`: Ultra-fast mode (step sizes 4 and 8 px, ~2.5s per image).
- `local_thresh_method` *(str)*: Local threshold algorithm: `"Li"`, `"Otsu"`, or `"Yen"`. Default: `"Li"`.
- `global_thresh_method` *(str)*: Global threshold algorithm: `"Otsu"`, `"Li"`, or `"Yen"`. Default: `"Otsu"`.
- `adjusting_factor_rng` *(tuple of floats)*: Dynamic factor range for adaptive sigmoid tuning curve. Default: `(-0.35, 0.25)`.

---

## 5. Understanding Output Files (4 Binary Masks)

For each processed image, the pipeline generates and saves exactly the **4 binary mask images** (saved with `dtype=bool`, matching `MitoSLIT_v4` format):

| Output Filename Suffix | Description |
| :--- | :--- |
| **`*(global).tif`** | Global threshold binary mask (without ROI masking). |
| **`*(global_TissueROI).tif`** | Global threshold binary mask restricted to the tissue ROI. |
| **`*(MitoSLIT_mask).tif`** | MitoSLIT locally iterative adaptive binary mask (without ROI masking). |
| **`*(MitoSLIT_mask_TissueROI).tif`** | MitoSLIT locally iterative adaptive binary mask restricted to the tissue ROI. |

*(If no tissue ROI is found, only the 2 unmasked images are saved).*

### Viewing Outputs in ImageJ / Fiji
- Because these are boolean masks (`True` / `False`), pixel values are `1` (mitochondria) and `0` (background).
- In Fiji/ImageJ, if the mask appears dark on initial open, simply press **Ctrl + Shift + C** (or go to *Image > Adjust > Brightness/Contrast*) and click **Auto** or set display range to `[0, 1]`.

---

## 6. Speed vs. Precision Modes

In `MitoSLITConfig`, you can choose between exact ground-truth precision and fast processing via `skipping_pixels_list`:

| Mode | `skipping_pixels_list` | Window Step Sizes | Runtime / Image | Output Consistency vs `v4` |
| :--- | :---: | :---: | :---: | :---: |
| **Standard Mode (Default)** | `(0, 1)` | 1 px, 2 px | ~35–120 s | **100% bit-for-bit identical** |
| **Fast Mode** | `(1, 3)` | 2 px, 4 px | ~9–15 s | **> 99.3% consistent** |
| **Ultra-Fast Mode** | `(3, 7)` | 4 px, 8 px | ~2–4 s | **> 99.1% consistent** |

To switch to Fast Mode, edit line 135 of [`MitoSLIT_pipeline.py`](./MitoSLIT_pipeline.py):
```python
skipping_pixels_list: Tuple[int, ...] = (1, 3)
```

---

## 7. Troubleshooting & FAQs

### Q: "No matching fluorescence images found"
- Verify that `input_dir` points to the correct folder.
- Check your `fluo_keyword`. For example, if your files are named `SAMPLE_Ch2_01.tif`, set `fluo_keyword = ("Ch2",)` or `fluo_keyword = "Ch2"`.
- Ensure your images are in `.tif` or `.tiff` format.

### Q: Too much background noise is segmented
- Increase `smallest_mito_size` (e.g. from `2` to `3`).
- Ensure `denoise_method = 1` is selected to subtract background autofluorescence.

### Q: Fine or faint mitochondrial branches are missing
- Decrease `tile_size` (e.g. try `10` or `15` for fine structures like neurons).
- Adjust the lower bound of `adjusting_factor_rng` (e.g. `(-0.45, 0.25)` to allow lower local thresholds in dim regions).
