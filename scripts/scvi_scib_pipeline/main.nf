// Pipeline to run scVI and SCIB with various parameters

process scVI {
    tag "Running scVI for sample ${meta.id} using parameters: batch_key=${batch_key}, n_hidden=${n_hidden}, n_latent=${n_latent}, n_layers=${n_layers}, dropout=${dropout}, lr=${lr}, batch_size=${batch_size}, likelihood=${likelihood}, n_epochs=${n_epochs}"
    container '/nfs/cellgeni/singularity/images/scvi-1.1.2.sif'
    input:
        tuple val(meta), path(adata)
        tuple val(batch_key), val(n_hidden), val(n_latent), val(n_layers), val(dropout), val(lr), val(batch_size), val(likelihood), val(n_epochs)

    output:
        tuple val(output_meta), path("*/representation.npy"), emit: npy
        tuple val(output_meta), path("*/model/model.pt"), emit: pt
        tuple val(output_meta), path("*/history.csv"), emit: csv
        tuple val(output_meta), path("*/summary.txt"), emit: txt
        path "versions.yml", emit: versions
    
    script:
        def args = task.ext.args ?: ""
        def param_string = "${batch_key}_${n_hidden}_${n_latent}_${n_layers}_${dropout}_${lr}_${batch_size}_${likelihood}_${n_epochs}"
        output_meta = [id: param_string, sample_id: meta.id]
        """
        # create output directory
        mkdir "$param_string"

        # run scVI integration
        integrate.py \
        $adata \
        "$param_string" \
        $args \
        --n_latent $n_latent \
        --n_layers $n_layers \
        --n_hidden $n_hidden \
        --learning_rate $lr \
        --batch_size $batch_size \
        --gene_likelihood $likelihood \
        --dropout_rate $dropout \
        --max_epochs $n_epochs \
        --batch_key $batch_key

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version)
            scvi: \$(python -c "import scvi; print(scvi.__version__)")
        END_VERSIONS
        """

}

process downsampleAdata {
    tag "Subsampling adata for sample ${meta.id}"
    container 'docker://quay.io/cellgeni/celltypist:1.6.3'

    input:
        tuple val(meta), path(adata)
    output:
        tuple val(meta), path("downsampled.h5ad"), emit: h5ad
        tuple val(meta), path("pca_matrix.npy"), emit: npy
        tuple val(meta), path("downsample_cells.json"), emit: json
        path "versions.yml", emit: versions
    
    script:
        def args = task.ext.args ?: ""
        """
        preprocess_downsample.py \
        $adata \
        $args \
        --output_indices downsample_cells.json \
        --output_adata downsampled.h5ad \
        --output_pca pca_matrix.npy

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version)
            scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
            celltypist: \$(python -c "import celltypist; print(celltypist.__version__)")
        END_VERSIONS
        """
}

// process downsample_adata {
//     tag "Subsampling adata for sample ${meta.id}"
//     container 'docker://quay.io/cellgeni/celltypist:1.6.3'

//     input:
//         tuple val(meta), path(adata)
//     output:
//         tuple val(meta), path("downsampled.h5ad"), emit: h5ad
//         tuple val(meta), path("downsample_cells.json"), emit: json
//         path "versions.yml", emit: versions
    
//     script:
//         def args = task.ext.args ?: ""
//         """
//         downsample_cells.py \
//         $adata \
//         $args \
//         --output_indices downsample_cells.json \
//         --output_adata downsampled.h5ad

//         cat <<-END_VERSIONS > versions.yml
//         "${task.process}":
//             python: \$(python --version)
//             celltypist: \$(python -c "import celltypist; print(celltypist.__version__)")
//         END_VERSIONS
//         """
// }

// process preprocess_pca {
//     tag "Preprocessing PCA for sample ${meta.id}"
//     container '/nfs/cellgeni/singularity/images/toh5ad.sif'

//     input:
//         tuple val(meta), path(adata)
//     output:
//         tuple val(output_meta), path("pca_matrix.npy"), emit: npy
//         path "versions.yml", emit: versions
    
//     script:
//         output_meta = [id: "PCA", sample_id: meta.id]
//         def args = task.ext.args ?: ""
//         """
//         preprocess_pca.py \
//         $adata \
//         $args \
//         --output_pca pca_matrix.npy

//         cat <<-END_VERSIONS > versions.yml
//         "${task.process}":
//             python: \$(python --version)
//             scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
//         END_VERSIONS
//         """
// }

process scIB {
    tag "Running scIB for sample ${meta_adata.id} and representation ${meta_representation.id}"
    container '/nfs/cellgeni/singularity/images/scib_metrics.sif'
    input:
        tuple val(meta_adata), path(adata)
        tuple val(meta_representation), path(representation)
        tuple val(meta_downsampled), path(downsampled_cells)
    output:
        tuple val(meta_representation), path("${meta_representation.id}.csv"), emit: csv
        path "versions.yml", emit: versions
    
    script:
        def args = task.ext.args ?: ""
        """
        benchmark.py \
        $adata \
        $representation \
        $downsampled_cells \
        $args \
        --representation_key ${meta_representation.id} \
        --n_jobs ${task.cpus} \
        --output "${meta_representation.id}.csv"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version)
            scib_metrics: \$(python -c "import scib_metrics; print(scib_metrics.__version__)")
        END_VERSIONS
        """
}

process scIB_plot {
    tag "Plotting SCIB results"
    container '/nfs/cellgeni/singularity/images/scib_metrics.sif'
    input:
        path result_files, stageAs: "results/*"
    output:
        path "*.svg", emit: svg
        path "versions.yml", emit: versions
    
    script:
        """
        plot_scib.py \
        --filedir results

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version)
            scib_metrics: \$(python -c "import scib_metrics; print(scib_metrics.__version__)")
        END_VERSIONS
        """
}

process collect_versions {
    tag "Collecting version files"
    input:
        path versions, stageAs: "versions/*.yml"
    output:
        path "versions.yml", emit: versions
    
    script:
        """
        # Collect versions from all processes
        cat versions/* > versions.yml
        """
}

workflow {
    sample_id = channel.value([[id: params.sample_id], file(params.adata)])

    // Convert the parameters to channels
    batch_key = channel.fromList(params.batch_key)
    n_hidden = channel.fromList(params.n_hidden)
    n_latent = channel.fromList(params.n_latent)
    n_layers = channel.fromList(params.n_layers)
    dropout = channel.fromList(params.dropout)
    lr = channel.fromList(params.lr)
    batch_size = channel.fromList(params.batch_size)
    likelihood = channel.fromList(params.likelihood)
    n_epochs = channel.fromList(params.n_epochs)

    // Combine all parameters into a single channel
    params = batch_key
                     .combine(n_hidden)
                     .combine(n_latent)
                     .combine(n_layers)
                     .combine(dropout)
                     .combine(lr)
                     .combine(batch_size)
                     .combine(likelihood)
                     .combine(n_epochs)

    // Run scVI
    scVI(sample_id, params)

    // Run downsampling
    downsampleAdata(sample_id)

    // Run PCA preprocessing
    // preprocess_pca(sample_id)

    // Run scIB
    downsampled_adata = downsampleAdata.output.h5ad.collect()
    downsampled_cells = downsampleAdata.output.json.collect()
    scIB(downsampled_adata, scVI.output.npy, downsampled_cells)

    // Plot scIB results
    scib_results = scIB.output.csv.map { meta, path -> path}.collect(flat: false)
    scIB_plot(scib_results)

    // Collect versions
    versions = scVI.output.versions.first().concat(downsampleAdata.output.versions, scIB.output.versions.first(), scIB_plot.output.versions).collect()
    collect_versions(versions)
}