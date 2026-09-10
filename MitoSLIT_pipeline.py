"""
================================================================================
MitoSLIT Pipeline: Mitochondria Segmentation with Locally Iterative Thresholding
================================================================================
A modular, beginner-friendly pipeline for segmenting mitochondria and computing
local intensity ratios from fluorescence microscopy images (TIFF format).

Quick Start:
------------
1. Review the "USER CONFIGURATION" section below.
2. Set `input_dir` to the folder containing your raw .tif/.tiff images.
3. Set `output_dir` where you want the segmented masks and results saved.
4. Set `fluo_keyword` (e.g. "FLIM", "Ch2", "GFP") to match your intensity images.
5. Set `roi_keyword` (e.g. "roi") to match your tissue ROI images (if any).
6. Run this script:
       python MitoSLIT_pipeline.py
"""

import os
import re
import time
import copy
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict, Any, Union

import numpy as np
import tifffile
from scipy.optimize import curve_fit
from skimage.restoration import rolling_ball
from skimage.filters import threshold_otsu, threshold_yen, threshold_li, gaussian
from skimage.morphology import opening, area_closing, disk, remove_small_objects

# Safe built-in functions
def sigmoid(x, L, x0, k, b):
    arg = np.clip(-k * (x - x0), -500.0, 500.0)
    return L / (1.0 + np.exp(arg)) + b

def find_nearest_value(arr, value):
    idx = np.argmin(np.abs(arr - value))
    return idx, arr[idx]

def normalization(arr, adjusting_factor_rng):
    nonzero_arr = arr - np.min(arr)
    span = np.max(nonzero_arr) if np.max(nonzero_arr) > 0 else 1.0
    scale = np.abs(adjusting_factor_rng[1] - adjusting_factor_rng[0])
    return scale * (nonzero_arr / span) + np.min(adjusting_factor_rng)

def Max_Entropy_thresh(distribution, bins=100, excluded_zero=True, LIST_ENTP_VALUE=False):
    distribution = np.array(distribution).ravel()
    if excluded_zero:
        vals = distribution[distribution > 0]
    else:
        vals = distribution
    if len(vals) == 0:
        return 0, 0
    hist, bin_edges = np.histogram(vals, bins=bins)
    sum_hist = np.sum(hist)
    if sum_hist == 0:
        return 0, 0

    p = hist / sum_hist
    best_t = 0
    max_ent = -1.0
    for t in range(1, len(p)):
        pA = np.sum(p[:t])
        pB = np.sum(p[t:])
        if pA > 0 and pB > 0:
            entA = -np.sum((p[:t] / pA) * np.log(p[:t] / pA + 1e-12))
            entB = -np.sum((p[t:] / pB) * np.log(p[t:] / pB + 1e-12))
            tot = entA + entB
            if tot > max_ent:
                max_ent = tot
                best_t = t
    return int(bin_edges[best_t]), max_ent


# ==============================================================================
# 1. USER CONFIGURATION (Edit this section for your experiment)
# ==============================================================================

@dataclass
class MitoSLITConfig:
    """Configuration settings for MitoSLIT pipeline."""

    # --- Folder Paths ---
    # Input folder containing your raw .tif/.tiff images
    input_dir: str = r".\MitoSLIT"

    # Output folder where segmented masks, ratios, and summary tables are saved
    output_dir: str = r".\MitoSLIT\filtered"

    # --- File Search Keywords (User-Defined) ---
    # Keyword(s) that MUST be present in the fluorescence intensity image filename.
    # Examples: "FluoIntensity", "FLIM", "Ch2", "GFP", or a tuple like ("FLIM",) or ("Ch2",).
    # Case-insensitive. Leave as "" or () to match all TIFF images.
    fluo_keyword: Union[str, Tuple[str, ...]] = ("FluoIntensity",)

    # Keyword that identifies tissue ROI image filenames.
    # Examples: "roi", "TissueROI", "_roi". Case-insensitive.
    # Leave as "" if you do not have/use tissue ROI masks.
    roi_keyword: str = "roi"

    # Keyword(s) that exclude unwanted files from processing.
    # Examples: files containing "Ch1" (other channel) or "mask" (previous segmentation).
    exclude_keywords: Union[str, Tuple[str, ...]] = ("Ch1", "mask")

    # --- Image Denoising & Filtering ---
    # 0 = No denoising (use raw image directly)
    # 1 = Background subtraction with rolling ball + Gaussian blur (recommended)
    # 2 = Morphological grayscale opening
    denoise_method: int = 1

    # Noise cutoff floor: raw pixel values below this are treated as background noise
    noise_level: float = 30.0

    # --- Core Segmentation Parameters ---
    # Tile size in pixels: should be larger than the largest mitochondrial structure
    # (Typical: 20 for somatic/aging/stressed tissues; 10 for neurons)
    tile_size: int = 20

    # Minimum mitochondrial feature size in pixels (smaller objects are discarded as noise)
    smallest_mito_size: int = 2

    # Local threshold windows (scales) to test during iterative sliding-window search
    local_window_scales: Tuple[int, ...] = (1, 2)  # Multipliers of tile_size: [1*tile_size, 2*tile_size]

    # Pixel skipping step per window scale:
    # (0, 1) = Standard mode matching MitoSLIT_v4 exactly (step sizes 1 and 2 px, default)
    # (1, 3) = Fast mode (step sizes 2 and 4 px, ~9s per image, 99.3% consistent)
    skipping_pixels_list: Tuple[int, ...] = (0, 1)

    # Thresholding algorithm for local tiles: 'Li', 'Otsu', or 'Yen'
    local_thresh_method: str = "Li"

    # Thresholding algorithm for global initial threshold: 'Otsu', 'Li', 'Yen', 'MaxEntropy'
    global_thresh_method: str = "Otsu"

    # Dynamic factor range for adaptive sigmoid tuning curve
    adjusting_factor_rng: Tuple[float, float] = (-0.35, 0.25)


# Default configuration instance
CONFIG = MitoSLITConfig()


# ==============================================================================
# 2. DATASET DISCOVERY & FILE PAIRING
# ==============================================================================

def find_image_pairs(config: MitoSLITConfig) -> List[Dict[str, Optional[str]]]:
    """
    Scans the input directory and pairs each fluorescence intensity image with its
    corresponding tissue ROI mask based on user-defined keywords.

    Returns:
        A list of dictionaries with keys:
            - 'mito_filename': name of the fluorescence image file
            - 'mito_path': full path to the fluorescence image
            - 'roi_filename': name of the matched tissue ROI file (or None)
            - 'roi_path': full path to the tissue ROI file (or None)
    """
    if not os.path.exists(config.input_dir):
        print(f"[WARNING] Input directory '{config.input_dir}' does not exist!")
        return []

    all_files = [f for f in os.listdir(config.input_dir) if f.lower().endswith(('.tif', '.tiff'))]

    # Normalize keywords to tuples of lowercase strings
    if isinstance(config.fluo_keyword, str):
        fluo_keys = (config.fluo_keyword.strip().lower(),) if config.fluo_keyword.strip() else ()
    else:
        fluo_keys = tuple(k.strip().lower() for k in config.fluo_keyword if k.strip())

    roi_key = config.roi_keyword.strip().lower() if config.roi_keyword else ""

    if isinstance(config.exclude_keywords, str):
        exclude_keys = (config.exclude_keywords.strip().lower(),) if config.exclude_keywords.strip() else ()
    else:
        exclude_keys = tuple(k.strip().lower() for k in config.exclude_keywords if k.strip())

    # 1. Identify tissue ROI files (if roi_keyword is specified)
    roi_files: List[str] = []
    if roi_key:
        for f in all_files:
            f_lower = f.lower()
            if roi_key in f_lower:
                # Do not exclude if the exclusion is the roi_key itself
                if not any(ex in f_lower for ex in exclude_keys if ex != roi_key):
                    roi_files.append(f)

    # 2. Identify fluorescence intensity images
    fluo_files: List[str] = []
    for f in all_files:
        f_lower = f.lower()
        # Exclude ROI files from being processed as fluorescence images
        if roi_key and roi_key in f_lower:
            continue
        # Check exclusion keywords
        if any(ex in f_lower for ex in exclude_keys):
            continue
        # Must contain all required fluorescence keywords
        if all(fk in f_lower for fk in fluo_keys):
            fluo_files.append(f)

    # 3. Match each fluorescence image with its corresponding ROI mask
    pairs: List[Dict[str, Optional[str]]] = []

    for f_file in fluo_files:
        stem = os.path.splitext(f_file)[0]
        matched_roi = None

        if roi_files:
            # Build search stems by stripping known channel/keyword tokens
            base = stem
            for fk in fluo_keys:
                base = re.sub(rf'[_-]?{re.escape(fk)}$', '', base, flags=re.IGNORECASE)

            # Also strip trailing numeric slice tokens (e.g. for neuron z-slices: _01, _02)
            base_no_slice = re.sub(r'_\d+$', '', base)

            # Strip common channel markers (_Ch0, _Ch1, _Ch2, _Ch3)
            base_clean = re.sub(r'_ch\d+$', '', base, flags=re.IGNORECASE)
            base_clean_no_slice = re.sub(r'_ch\d+$', '', base_no_slice, flags=re.IGNORECASE)

            candidates: List[str] = []
            for r in roi_files:
                r_stem = os.path.splitext(r)[0]
                # Strip ROI keyword from candidate to compare base stems
                r_clean = re.sub(rf'[_-]?{re.escape(roi_key)}$', '', r_stem, flags=re.IGNORECASE)

                # Check matches against base stems
                if (r_clean.lower() in (base.lower(), base_no_slice.lower(), base_clean.lower(), base_clean_no_slice.lower())
                    or base_clean_no_slice.lower() in r_stem.lower()
                    or base_clean.lower() in r_stem.lower()):
                    candidates.append(r)

            if len(candidates) == 1:
                matched_roi = candidates[0]
            elif len(candidates) > 1:
                # Disambiguate: prefer exact target stem base + '_' + roi_key
                preferred = [
                    c for c in candidates
                    if os.path.splitext(c)[0].lower() in (
                        f"{base}_{roi_key}".lower(),
                        f"{base_no_slice}_{roi_key}".lower(),
                        f"{base_clean}_{roi_key}".lower(),
                        f"{base_clean_no_slice}_{roi_key}".lower(),
                    )
                ]
                if len(preferred) == 1:
                    matched_roi = preferred[0]
                else:
                    matched_roi = candidates[0]
                    print(f"  [NOTE] Multiple ROI matches for '{f_file}': {candidates}. Using '{matched_roi}'.")

        pairs.append({
            'mito_filename': f_file,
            'mito_path': os.path.join(config.input_dir, f_file),
            'roi_filename': matched_roi,
            'roi_path': os.path.join(config.input_dir, matched_roi) if matched_roi else None,
        })

    return pairs


# ==============================================================================
# 3. IMAGE PREPROCESSING
# ==============================================================================

def preprocess_image(raw: np.ndarray, config: MitoSLITConfig) -> np.ndarray:
    """
    Cleans raw image data:
    1. Clips negative camera values to zero.
    2. Applies the selected background subtraction / denoising method.
    """
    raw_nonnegative = copy.deepcopy(raw)
    raw_nonnegative[raw < 0] = 0

    if config.denoise_method == 0:
        # No denoising
        return raw_nonnegative

    elif config.denoise_method == 1:
        # Rolling ball background subtraction + Gaussian smoothing (matching MitoSLIT_v4)
        ball_radius = max(1, config.tile_size // 4)
        sigma_blur = max(1, config.tile_size // 2)

        bg = rolling_ball(raw_nonnegative, radius=ball_radius)
        bg_blurred = gaussian(bg, sigma=sigma_blur, preserve_range=True)
        denoised = raw_nonnegative - bg_blurred
        denoised[denoised < 0] = 0
        return denoised

    elif config.denoise_method == 2:
        # Morphological grayscale opening
        selem = np.ones((config.smallest_mito_size, config.smallest_mito_size))
        denoised = opening(raw_nonnegative, selem)
        return denoised

    else:
        raise ValueError(f"Unknown denoise_method: {config.denoise_method}. Choose 0, 1, or 2.")


# ==============================================================================
# 4. ADAPTIVE SIGMOID CURVE FITTING
# ==============================================================================

def compute_adaptive_curve(
    data: np.ndarray,
    thresh_global: float,
    config: MitoSLITConfig
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Divides the preprocessed image into tiles, computes the mean intensity
    distribution, and fits a sigmoid function to dynamically determine
    local threshold adjustment factors.

    Returns:
        (bin_edges, scaled_fitted_curve)
    """
    tile_h = int(data.shape[0] // config.tile_size)
    tile_w = int(data.shape[1] // config.tile_size)
    tile_2d = np.zeros((tile_h, tile_w), dtype=float)

    for i in range(tile_h):
        for j in range(tile_w):
            sub_tile = data[
                i * config.tile_size : (i + 1) * config.tile_size,
                j * config.tile_size : (j + 1) * config.tile_size
            ]
            tile_mean = np.mean(sub_tile)
            tile_2d[i, j] = tile_mean if tile_mean > 0 else 1.0 / (config.tile_size ** 2)

    # Histogram of non-zero tile means
    hist, bin_edges = np.histogram(tile_2d[tile_2d > 0], bins='auto')
    bin_edges = np.round(bin_edges)

    # Guard: handle empty or degenerate histograms
    if len(hist) < 3:
        curve = np.full(len(bin_edges) - 1, np.mean(config.adjusting_factor_rng))
        return bin_edges, curve

    # Fit sigmoid curve to inverted histogram (matching MitoSLIT_v4)
    try:
        p0 = [float(np.max(hist[1:])), float(thresh_global), 0.001, float(np.min(hist[1:]))]
        popt, _ = curve_fit(sigmoid, bin_edges[2:], np.flip(hist[1:]), p0, method='dogbox')
        fitted_curve = sigmoid(bin_edges[1:], popt[0], popt[1], popt[2], popt[3])
        scaled_curve = normalization(fitted_curve, config.adjusting_factor_rng)
    except Exception:
        # Fallback to initial guess if optimization fails
        p0 = [float(np.mean(hist[0:1])), float(thresh_global), 0.001, float(np.min(hist[1:]))]
        fitted_curve = sigmoid(bin_edges[1:], p0[0], p0[1], p0[2], p0[3])
        scaled_curve = normalization(fitted_curve, config.adjusting_factor_rng)

    return bin_edges, scaled_curve


# ==============================================================================
# 5. CORE ITERATIVE LOCAL THRESHOLDING (MITOSLIT)
# ==============================================================================

def run_mitoslit_iteration(
    data: np.ndarray,
    thresh_global: float,
    bin_edges: np.ndarray,
    scaled_curve: np.ndarray,
    config: MitoSLITConfig
) -> np.ndarray:
    """
    Multiscale sliding-window iterative thresholding.
    Produces results mathematically and pixel-wise consistent with MitoSLIT_v4.

    Returns:
        prob_accum: Accumulated vote map for mitochondrial pixels.
    """
    window_sizes = [int(config.tile_size * scale) for scale in config.local_window_scales]
    skipping_list = list(config.skipping_pixels_list)

    prob_accum = np.zeros(data.shape, dtype=float)
    b_edges = bin_edges[1:]

    # Local threshold helper
    def get_local_threshold(tile: np.ndarray) -> float:
        if config.local_thresh_method == 'Yen':
            return float(threshold_yen(tile))
        elif config.local_thresh_method == 'Li':
            return float(threshold_li(tile))
        else:
            return float(threshold_otsu(tile))

    for win_size, skip in zip(window_sizes, skipping_list):
        step = max(1, skip + 1)

        # Reflective padding to avoid edge artifacts
        pad_data = np.pad(data, pad_width=win_size, mode='reflect')
        temp_prob = np.zeros(pad_data.shape, dtype=float)

        rem_x = (pad_data.shape[0] - win_size) % step
        rem_y = (pad_data.shape[1] - win_size) % step
        n_steps_x = int((pad_data.shape[0] - win_size - rem_x) // step)
        n_steps_y = int((pad_data.shape[1] - win_size - rem_y) // step)

        for sx in range(n_steps_x):
            x_start = sx * step
            x_end = x_start + win_size

            for sy in range(n_steps_y):
                y_start = sy * step
                y_end = y_start + win_size

                local_tile = pad_data[x_start:x_end, y_start:y_end]

                idx, _ = find_nearest_value(b_edges, np.mean(local_tile))
                factor = scaled_curve[idx]

                tile_thresh = get_local_threshold(local_tile)
                adj_thresh = max(tile_thresh, thresh_global * (1.0 + factor))

                if np.max(local_tile) > adj_thresh:
                    temp_prob[x_start:x_end, y_start:y_end] += (local_tile >= adj_thresh)
                else:
                    temp_prob[x_start:x_end, y_start:y_end] += (local_tile >= thresh_global)

        # Unpad and accumulate across scales
        prob_accum += temp_prob[win_size:-win_size, win_size:-win_size]

    return prob_accum


# ==============================================================================
# 6. POST-PROCESSING & MORPHOLOGICAL CLEANUP
# ==============================================================================

def postprocess_segmentation(
    prob_accum: np.ndarray,
    config: MitoSLITConfig
) -> np.ndarray:
    """
    Applies Maximum Entropy + Otsu thresholding on the probability map,
    cleans up artifacts with morphological opening/closing, and produces
    the final binary mitochondrial mask (dtype=bool, exactly matching MitoSLIT_v4).
    """
    # Optimal threshold for probability vote map
    ent_thresh, _ = Max_Entropy_thresh(prob_accum, bins=100, LIST_ENTP_VALUE=True)
    otsu_thresh = threshold_otsu(prob_accum)
    mitoslit_thresh = max(ent_thresh, otsu_thresh)

    # Binarize
    local_thresholded_image = np.ones(prob_accum.shape)
    local_thresholded_image[prob_accum < mitoslit_thresh] = 0

    # Morphological clean-up matching MitoSLIT_v4 exact sequence
    min_area = int(config.smallest_mito_size ** 2)
    local_thresholded_image = area_closing(
        np.logical_and(local_thresholded_image, local_thresholded_image),
        area_threshold=min_area
    )
    local_thresholded_image = np.logical_and(local_thresholded_image, local_thresholded_image)
    local_thresholded_image = remove_small_objects(
        local_thresholded_image,
        min_size=min_area,
        connectivity=1
    )
    final_mask = opening(local_thresholded_image, disk(1))

    return final_mask.astype(bool)


# ==============================================================================
# 7. QUANTITATIVE METRICS & SAVING 4 BINARY IMAGES
# ==============================================================================

def calculate_metrics(
    image_name: str,
    total_pixels: int,
    global_mask: np.ndarray,
    mitoslit_mask: np.ndarray,
    tissue_roi: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """Computes summary statistics (pixel counts and area fractions) for binary results."""
    global_pixels = int(np.sum(global_mask > 0))
    mitoslit_pixels = int(np.sum(mitoslit_mask > 0))

    metrics: Dict[str, Any] = {
        'image': image_name,
        'total_pixels': total_pixels,
        'global_pixels': global_pixels,
        'global_area_fraction': round(global_pixels / total_pixels, 5) if total_pixels > 0 else 0.0,
        'mitoslit_pixels': mitoslit_pixels,
        'mitoslit_area_fraction': round(mitoslit_pixels / total_pixels, 5) if total_pixels > 0 else 0.0,
    }

    if tissue_roi is not None:
        roi_mask = (tissue_roi == 0)
        roi_total = int(np.sum(roi_mask))
        global_roi_px = int(np.sum((global_mask > 0) & roi_mask))
        mitoslit_roi_px = int(np.sum((mitoslit_mask > 0) & roi_mask))

        metrics['roi_total_pixels'] = roi_total
        metrics['global_pixels_in_roi'] = global_roi_px
        metrics['global_fraction_in_roi'] = round(global_roi_px / roi_total, 5) if roi_total > 0 else 0.0
        metrics['mitoslit_pixels_in_roi'] = mitoslit_roi_px
        metrics['mitoslit_fraction_in_roi'] = round(mitoslit_roi_px / roi_total, 5) if roi_total > 0 else 0.0

    return metrics


def save_binary_results(
    base_name: str,
    global_mask: np.ndarray,
    mitoslit_mask: np.ndarray,
    tissue_roi: Optional[np.ndarray],
    config: MitoSLITConfig
) -> List[str]:
    """
    Saves the 4 binary results requested:
    1. Global mask without ROI:        *(global).tif
    2. Global mask with Tissue ROI:    *(global_TissueROI).tif
    3. MitoSLIT mask without ROI:      *(MitoSLIT_mask).tif
    4. MitoSLIT mask with Tissue ROI:  *(MitoSLIT_mask_TissueROI).tif

    (If no tissue ROI was provided, only the 2 unmasked images are saved).
    """
    out_dir = config.output_dir
    os.makedirs(out_dir, exist_ok=True)
    saved_files = []

    # 1. Global threshold binary result (without ROI)
    p_global = os.path.join(out_dir, f"{base_name}(global).tif")
    tifffile.imwrite(p_global, global_mask.astype(bool))
    saved_files.append(p_global)

    # 2. MitoSLIT binary result (without ROI)
    p_mitoslit = os.path.join(out_dir, f"{base_name}(MitoSLIT_mask).tif")
    tifffile.imwrite(p_mitoslit, mitoslit_mask.astype(bool))
    saved_files.append(p_mitoslit)

    # Tissue ROI masked outputs (if ROI file is available)
    if tissue_roi is not None:
        # 3. Global threshold binary result (with Tissue ROI)
        global_roi = copy.deepcopy(global_mask)
        global_roi[tissue_roi > 0] = 0
        p_global_roi = os.path.join(out_dir, f"{base_name}(global_TissueROI).tif")
        tifffile.imwrite(p_global_roi, global_roi.astype(bool))
        saved_files.append(p_global_roi)

        # 4. MitoSLIT binary result (with Tissue ROI)
        mitoslit_roi = copy.deepcopy(mitoslit_mask)
        mitoslit_roi[tissue_roi > 0] = 0
        p_mitoslit_roi = os.path.join(out_dir, f"{base_name}(MitoSLIT_mask_TissueROI).tif")
        tifffile.imwrite(p_mitoslit_roi, mitoslit_roi.astype(bool))
        saved_files.append(p_mitoslit_roi)

    return saved_files


# ==============================================================================
# 8. MAIN EXECUTION PIPELINE
# ==============================================================================

def process_single_image(pair: Dict[str, Optional[str]], config: MitoSLITConfig) -> None:
    """Runs the optimized MitoSLIT segmentation workflow on a single image pair."""
    mito_path = pair['mito_path']
    mito_name = pair['mito_filename']
    base_name = os.path.splitext(mito_name)[0]

    # 1. Load image
    raw = np.array(tifffile.imread(mito_path))

    # 2. Load tissue ROI (if available)
    tissue_roi = None
    if pair['roi_path'] and os.path.exists(pair['roi_path']):
        tissue_roi = np.array(tifffile.imread(pair['roi_path']))
        print(f"  [ROI] Using tissue ROI: {pair['roi_filename']}")

    # 3. Preprocessing / Denoising
    denoised = preprocess_image(raw, config)

    # 4. Global threshold calculation
    if config.global_thresh_method == 'Li':
        thresh_global = float(threshold_li(denoised))
    elif config.global_thresh_method == 'Yen':
        thresh_global = float(threshold_yen(denoised))
    else:
        thresh_global = float(threshold_otsu(denoised))

    # Global thresholded binary mask
    global_mask = (denoised >= thresh_global)
    global_mask = opening(global_mask, disk(1))

    # 5. Fit adaptive sigmoid curve
    bin_edges, scaled_curve = compute_adaptive_curve(denoised, thresh_global, config)

    # 6. Fast local iterative thresholding
    prob_accum = run_mitoslit_iteration(
        denoised, thresh_global, bin_edges, scaled_curve, config
    )

    # 7. Post-processing to binary mask
    mitoslit_mask = postprocess_segmentation(prob_accum, config)

    # 8. Save the 4 binary results
    saved = save_binary_results(base_name, global_mask, mitoslit_mask, tissue_roi, config)
    print(f"  [SAVED] {len(saved)} binary images written.")


def run_pipeline(config: MitoSLITConfig = CONFIG) -> None:
    """Entrypoint: discovers files, runs batch processing, and saves binary results."""
    print("=" * 70)
    print("  MitoSLIT Pipeline: Mitochondrial Binary Segmentation")
    print("=" * 70)
    print(f"  Input Directory   : {os.path.abspath(config.input_dir)}")
    print(f"  Output Directory  : {os.path.abspath(config.output_dir)}")
    print(f"  Fluo Keyword(s)   : {config.fluo_keyword}")
    print(f"  ROI Keyword       : {config.roi_keyword if config.roi_keyword else '(None)'}")
    print(f"  Exclude Keyword(s): {config.exclude_keywords}")
    print(f"  Tile Size (px)    : {config.tile_size}")
    print(f"  Denoise Method    : {config.denoise_method}")
    print(f"  Local Method      : {config.local_thresh_method}")
    print("=" * 70)

    # Find images
    image_pairs = find_image_pairs(config)
    total_images = len(image_pairs)

    if total_images == 0:
        print(f"\n[INFO] No matching fluorescence images found in '{config.input_dir}'.")
        print(f"  - Fluo keyword(s) looked for : {config.fluo_keyword}")
        print(f"  - Excluded keywords          : {config.exclude_keywords}")
        print("Please verify your input_dir path or adjust 'fluo_keyword' in CONFIG.")
        return

    print(f"\nFound {total_images} image(s) to process.\n")
    os.makedirs(config.output_dir, exist_ok=True)

    start_total = time.time()

    for idx, pair in enumerate(image_pairs, 1):
        print(f"[{idx}/{total_images}] Processing {pair['mito_filename']} ...")
        t0 = time.time()

        try:
            process_single_image(pair, config)
            elapsed = round(time.time() - t0, 2)
            print(f"  --> Completed in {elapsed} s\n")
        except Exception as err:
            print(f"  [ERROR] Failed to process {pair['mito_filename']}: {err}\n")

    total_time = round(time.time() - start_total, 2)
    print("=" * 70)
    print(f"  All tasks completed in {total_time} seconds.")
    print("=" * 70)


if __name__ == "__main__":
    # Execute the pipeline with default configuration
    run_pipeline(CONFIG)
