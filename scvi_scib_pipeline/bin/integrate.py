#!/usr/bin/env python3
import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import scvi

SCVI_LATENT_KEY = "X_scVI"


def init_parser() -> argparse.ArgumentParser:
    """
    Initialise argument parser for the script
    """
    parser = argparse.ArgumentParser(
        description="A script for data integration using scVI"
    )
    parser.add_argument(
        "input",
        metavar="file",
        type=str,
        help="specify a path to input file",
    )
    parser.add_argument(
        "output",
        metavar="file",
        type=str,
        help="specify a path to output file",
        default="results",
    )
    parser.add_argument(
        "--learning_rate",
        metavar="float",
        type=float,
        help="specify the learning rate",
        default=0.0001,
    )
    parser.add_argument(
        "--batch_size",
        metavar="int",
        type=int,
        help="specify the batch size",
        default=128,
    )
    parser.add_argument(
        "--n_latent",
        metavar="int",
        type=int,
        help="specify the number of latent dimensions",
        default=10,
    )
    parser.add_argument(
        "--n_layers",
        metavar="int",
        type=int,
        help="specify the number of layers",
        default=2,
    )
    parser.add_argument(
        "--n_hidden",
        metavar="int",
        type=int,
        help="specify the number of hidden units",
        default=128,
    )
    parser.add_argument(
        "--gene_likelihood",
        metavar="str",
        type=str,
        help="specify the gene likelihood",
        default="nb",
    )
    parser.add_argument(
        "--dropout_rate",
        metavar="float",
        type=float,
        help="specify the dropout rate",
        default=0.1,
    )
    parser.add_argument(
        "--max_epochs",
        metavar="int",
        type=int,
        help="specify the maximum number of epochs",
        default=None,
    )
    parser.add_argument(
        "--batch_key",
        metavar="str",
        type=str,
        help="specify batch key for integration",
        default="donor_batch",
    )
    parser.add_argument(
        "--layer",
        metavar="str",
        type=str,
        help="specify the layer to use for training",
        default=None,
    )
    return parser


def train_model(
    adata,
    layer=None,
    batch_key="donor_batch",
    n_latent=30,
    n_layers=2,
    n_hidden=128,
    dropout_rate=0.1,
    gene_likelihood="nb",
    max_epochs=None,
    batch_size=128,
    lr=0.0001,
):
    # Set the seed
    scvi.settings.seed = 4

    # Set up the anndata object
    scvi.model.SCVI.setup_anndata(adata, layer=layer, batch_key=batch_key)

    # Create a model object
    model = scvi.model.SCVI(
        adata=adata,
        n_layers=n_layers,
        n_latent=n_latent,
        n_hidden=n_hidden,
        dropout_rate=dropout_rate,
        gene_likelihood=gene_likelihood,
    )

    # Train the model
    plan_kwargs = {"lr": lr}
    model.train(max_epochs=max_epochs, batch_size=batch_size, plan_kwargs=plan_kwargs)

    return model


def save_model(
    model,
    output_dir="results",
):
    os.makedirs(output_dir, exist_ok=True)

    # Save representation
    print("Saving representation")
    representation = model.get_latent_representation()
    np.save(
        os.path.join(output_dir, "representation.npy"),
        representation,
    )

    # Save model
    print("Saving model")
    if os.path.exists(os.path.join(output_dir, "model")):
        os.remove(os.path.join(output_dir, "model", "model.pt"))
        os.rmdir(os.path.join(output_dir, "model"))

    model.save(os.path.join(output_dir, "model"))

    # Save summary
    print("Saving summary")
    with open(os.path.join(output_dir, "summary.txt"), "w") as file:
        file.write(model.summary_string)
        file.write(
            f"reconstruction_error={model.get_reconstruction_error(adata)['reconstruction_loss']}"
        )

    # Save training history
    print("Saving training history")
    history = pd.concat(model.history, axis=1)
    history.columns = history.columns.droplevel(level=0)
    history.to_csv(os.path.join(output_dir, "history.csv"))


if __name__ == "__main__":
    # parse script arguments
    parser = init_parser()
    args = parser.parse_args()

    print(f"Running integration with following params: {args}")

    # read data
    adata = sc.read(args.input)

    # train model
    model = train_model(
        adata,
        layer=args.layer,
        batch_key=args.batch_key,
        n_latent=args.n_latent,
        n_layers=args.n_layers,
        n_hidden=args.n_hidden,
        dropout_rate=args.dropout_rate,
        gene_likelihood=args.gene_likelihood,
        max_epochs=args.max_epochs,
        batch_size=args.batch_size,
        lr=args.learning_rate,
    )

    # save model
    save_model(
        model,
        output_dir=args.output,
    )
