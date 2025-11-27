process scIBenchmark {
    tag "Running scIB for sample ${meta_adata.id} and representation ${meta_representation.id}"
    container 'docker://quay.io/cellgeni/scib-metrics:latest'
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
            jax: \$(python -c "import jax; print(jax.__version__)")
        END_VERSIONS
        """
}

process scIBPlot {
    tag "Plotting SCIB results"
    container '/nfs/cellgeni/singularity/images/scib_metrics.sif'
    input:
        path result_files, stageAs: "results/*"
    output:
        path "*.svg", emit: svg
        path "*.csv", emit: csv
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

workflow {
    representations = Channel.fromPath(params.representations, checkIfExists: true)
        .splitCsv(header: true, sep: ',')
        .map { row -> tuple(row, file(row.path)) }
    
    downsample_adata = Channel.value([[id: "scib_run"], file(params.downsample_adata)])
    downsample_cells = Channel.value([[id: "scib_run"], file(params.downsample_cells)])

    // representations.take(5).view()
    // downsample_adata.view()
    // downsample_cells.view()

    scIBenchmark(
        downsample_adata,
        representations,
        downsample_cells
    )

    scib_results = scIBenchmark.output.csv.map { _meta, path -> path}.collect(flat: false)
    scIBPlot(scib_results)


}