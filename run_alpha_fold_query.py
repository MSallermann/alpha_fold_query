from pathlib import Path
import pandas as pd
import alpha_fold_query
import logging
from collections.abc import Sequence


logger = logging.getLogger(__name__)


def main(uniprot_ids_in: Sequence[str], ouput_path: Path):

    if output_path.exists():
        df_out = pd.read_parquet(output_path, engine="pyarrow")
    else:
        df_out = pd.DataFrame(
            columns=["uniprot_id", "n_residues", "sequence", "plddts"]
        )

    uniprot_ids_out_list = list(df_out["uniprot_id"])
    seen = set(uniprot_ids_out_list)

    # Remove all the ids we have seen already
    uniprot_ids_in_unique = set(uniprot_ids_in)
    ids_to_query = uniprot_ids_in_unique.difference(seen)

    logger.info(
        f"Of {len(uniprot_ids_in)} input ids, {len(uniprot_ids_in_unique)} are unique. Out of these {len(seen)} are already in the output file '{ouput_path}'."
    )
    logger.info(f"Therefore I will query {len(ids_to_query)} ids.")

    query_result_generator = alpha_fold_query.query_alphafold_bulk(
        list(ids_to_query), retries=1, backoff_time=1
    )

    idx = 0
    try:
        for query_result in query_result_generator:
            idx += 1

            logger.info(
                f"Queried {query_result.accession} [{idx} / {len(ids_to_query)}]"
            )

            if query_result.http_status != 200:
                continue

            df_out.loc[len(df_out)] = [
                query_result.accession,
                len(query_result.sequence),
                query_result.sequence,
                query_result.plddts,
            ]

    except BaseException as e:
        output_path_exc = output_path.with_name("saved_after_exc").with_suffix(
            ".parquet"
        )

        logger.exception(
            f"Encountered exception {e}. Will try to save data to {output_path_exc}",
            stack_info=True,
            stacklevel=1,
        )

        df_out.to_parquet(output_path_exc, index=False, engine="pyarrow")

        raise e

    logger.info(f"Saving to {output_path}")
    df_out.to_parquet(output_path)


if __name__ == "__main__":
    from rich.logging import RichHandler

    FORMAT = "%(message)s"
    logging.basicConfig(
        level=logging.INFO, format=FORMAT, datefmt="[%X]", handlers=[RichHandler()]
    )
    data_path = Path("./IDRome.csv")
    output_path = Path("./data_out.parquet")

    df_in = pd.read_csv(data_path)

    ids_to_query = list(df_in["UniProt_ID"])

    main(uniprot_ids_in=ids_to_query, ouput_path=output_path)
