# Split and Combined Downsampling Pipeline

This directory contains multiple scripts for preprocessing, downsampling, and benchmarking single-cell data:

## Scripts

### 1. preprocess_downsample.py (Combined Script)

**NEW: Combined script that performs PCA calculation followed by downsampling in a single run.**

**Functionality:**
- Loads AnnData object
- Optional normalization and log transformation
- Optional highly variable gene detection
- PCA calculation
- Cell downsampling based on cell types
- Saves PCA matrix as `.npy` file
- Saves downsampled cell indices in multiple formats

**Usage:**
```bash
./preprocess_downsample.py input_file.h5ad --celltype_key CELLTYPE --fraction FRACTION [options]
```

**Key Arguments:**
- `input_file`: Path to input AnnData file (.h5ad)
- `--celltype_key`: Column name in adata.obs containing cell type annotations (required)
- `--fraction`: Fraction to downsample to, between 0 and 1 (required)
- `--lognorm`: Apply normalization and log1p transformation
- `--hvg`: Find highly variable genes before PCA
- `--n_hvg`: Number of highly variable genes (default: 2000)
- `--n_pca`: Number of PCA components (default: 50)
- `--output_pca`: Output path for PCA matrix (default: pca.npy)
- `--output_indices`: Output file for cell indices (default: downsampled_indices.json)
- `--indices_format`: Format for indices (json, txt, csv) (default: json)
- `--output_adata`: Optional path to save processed AnnData
- `--log_level`: Logging level

**Example:**
```bash
./preprocess_downsample.py data.h5ad --celltype_key cell_type --fraction 0.1 --lognorm --hvg --output_pca my_pca.npy --output_indices my_indices.json
```

### 2. preprocess_pca.py (Individual Script)

Performs data preprocessing and generates PCA matrix.

**Functionality:**
- Loads AnnData object
- Optional normalization and log transformation
- Optional highly variable gene detection
- PCA calculation
- Saves PCA matrix as `.npy` file
- Optionally saves processed AnnData object

**Usage:**
```bash
./preprocess_pca.py input_file.h5ad [options]
```

**Key Arguments:**
- `input_file`: Path to input AnnData file (.h5ad)
- `--lognorm`: Apply normalization and log1p transformation
- `--hvg`: Find highly variable genes before PCA
- `--n_hvg`: Number of highly variable genes (default: 2000)
- `--n_pca`: Number of PCA components (default: 50)
- `--output_pca`: Output path for PCA matrix (default: pca.npy)
- `--output_adata`: Optional path to save processed AnnData
- `--log_level`: Logging level (DEBUG, INFO, WARNING, ERROR)

**Output Files:**
- `pca.npy`: PCA matrix as numpy array
- `pca_barcodes.txt`: Cell barcodes corresponding to PCA rows
- `preprocess_pca.log`: Log file

**Example:**
```bash
./preprocess_pca.py data.h5ad --lognorm --hvg --n_hvg 3000 --n_pca 100 --output_pca my_pca.npy
```

### 3. downsample_cells.py (Individual Script)

Performs cell downsampling and outputs selected cell indices.

**Functionality:**
- Loads AnnData object
- Performs stratified downsampling by cell type
- Returns indices of selected cells
- Multiple output formats supported (JSON, TXT, CSV)
- Optionally saves downsampled AnnData object

**Usage:**
```bash
./downsample_cells.py input_file.h5ad --celltype_key CELLTYPE --fraction 0.1 [options]
```

**Key Arguments:**
- `input_file`: Path to input AnnData file (.h5ad)
- `--celltype_key`: Column name in adata.obs containing cell type annotations (required)
- `--fraction`: Fraction to downsample to, between 0 and 1 (required)
- `--output_indices`: Output file for cell indices (default: downsampled_indices.json)
- `--indices_format`: Format for indices (json, txt, csv) (default: json)
- `--output_adata`: Optional path to save downsampled AnnData
- `--log_level`: Logging level (DEBUG, INFO, WARNING, ERROR)

**Output Formats:**
1. **JSON** (default): Complete information with metadata
   ```json
   {
     "selected_indices": [0, 15, 23, ...],
     "selected_barcodes": ["CELL_0001", "CELL_0015", ...],
     "total_selected": 1000,
     "metadata": {
       "description": "Downsampled cell indices and barcodes",
       "format": "JSON"
     }
   }
   ```

2. **TXT**: Simple list of cell barcodes (one per line)
   ```
   CELL_0001
   CELL_0015
   CELL_0023
   ...
   ```

3. **CSV**: Table with original indices and barcodes
   ```csv
   original_index,barcode
   0,CELL_0001
   15,CELL_0015
   23,CELL_0023
   ...
   ```

**Output Files:**
- Indices file (format specified by `--indices_format`)
- `downsample_cells.log`: Log file

**Example:**
```bash
./downsample_cells.py data.h5ad --celltype_key cell_type --fraction 0.1 --output_indices selected_cells.json --indices_format json
```

### 4. benchmark.py

Runs scib benchmark on downsampled AnnData with representation matrix and pre-integrated PCA.

**Functionality:**
- Loads downsampled AnnData object (with PCA), representation matrix, and downsampled indices
- Attaches representation matrix to the downsampled AnnData using selected cell indices
- Runs scib benchmark with bio-conservation and batch-correction metrics
- Uses X_pca as pre-integrated embedding for baseline comparison
- Saves benchmark results

**Usage:**
```bash
./benchmark.py downsampled_adata.h5ad representation.npy downsampled_indices.json --batch_key BATCH --celltype_key CELLTYPE [options]
```

**Key Arguments:**
- `downsampled_adata`: Path to downsampled AnnData file (.h5ad) with PCA
- `representation_file`: Path to full representation matrix (.npy) from original dataset
- `downsampled_file`: Path to downsampled indices (.json)
- `--batch_key`: Column name in adata.obs containing batch information (required)
- `--celltype_key`: Column name in adata.obs containing cell type annotations (required)
- `--representation_key`: Key to store representation in adata.obsm (default: X_benchmark)
- `--embedding_keys`: List of embedding keys to benchmark (if not provided, uses representation_key)
- `--n_jobs`: Number of parallel jobs (default: 1)
- `--output`: Output CSV file for results (default: benchmark_results.csv)
- `--log_level`: Logging level

**Output Files:**
- `benchmark_results.csv`: Benchmark metrics results
- `benchmark.log`: Log file

**Example:**
```bash
./benchmark.py downsampled_data.h5ad pca_matrix.npy selected_cells.json --batch_key batch --celltype_key cell_type --n_jobs 4 --output my_benchmark.csv
```

## Pipeline Usage Options

### Option 1: Combined Script (Recommended for most use cases)

Run preprocessing, PCA calculation, and downsampling in one step:

```bash
./preprocess_downsample.py data.h5ad --celltype_key cell_type --fraction 0.1 --lognorm --hvg --output_pca pca.npy --output_indices indices.json
```

Then run benchmarking:

```bash
./benchmark.py downsampled_output.h5ad pca.npy indices.json --batch_key batch --celltype_key cell_type --output benchmark_results.csv
```

### Option 2: Individual Scripts (For maximum flexibility)

1. **Preprocessing step:**
   ```bash
   ./preprocess_pca.py original_data.h5ad --lognorm --hvg --output_pca pca_matrix.npy --output_adata processed_data.h5ad
   ```

2. **Downsampling step:**
   ```bash
   ./downsample_cells.py processed_data.h5ad --celltype_key cell_type --fraction 0.1 --output_indices downsampled_indices.json
   ```

3. **Benchmarking step:**
   ```bash
   ./benchmark.py processed_downsampled_data.h5ad pca_matrix.npy downsampled_indices.json --batch_key batch --celltype_key cell_type --output benchmark_results.csv
   ```

## Dependencies

All scripts require:
- `anndata`
- `scanpy`
- `celltypist`
- `numpy`

Additional dependencies for benchmarking:
- `pandas`
- `scib-metrics`

## Log Files

Each script generates its own log file:
- `preprocess_downsample.log` (combined script)
- `preprocess_pca.log` (individual preprocessing)
- `downsample_cells.log` (individual downsampling)
- `benchmark.log` (benchmarking)

Use `--log_level DEBUG` for more detailed logging.

## Notes

- The preprocessing script saves cell barcodes alongside the PCA matrix for reference
- The downsampling script uses `celltypist.samples.downsample_adata` with stratified sampling
- The benchmark script automatically subsets both the AnnData and representation matrix to match the downsampled cells
- Random seed is set to 4 for reproducible downsampling
- All scripts include comprehensive logging and error handling
- The benchmark script uses scib-metrics with both bio-conservation and batch-correction metrics
