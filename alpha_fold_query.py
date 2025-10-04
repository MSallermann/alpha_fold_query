from typing import Generator, Any
import requests
import time
from collections.abc import Sequence
from collections import namedtuple
import logging

logger = logging.getLogger()

# AlphaFold DB API: returns JSON metadata with model file URLs keyed by UniProt accession
# e.g. https://alphafold.ebi.ac.uk/api/prediction/P69905
ALPHAFOLD_PREDICTION_URL = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"


def get_with_retries(
    url, retries: int, backoff_time: int, **kwargs
) -> requests.Response:

    for i in range(retries):
        r = requests.get(url, **kwargs)
        if r.status_code != 200:
            time.sleep(backoff_time)
            continue
        else:
            return r

    return requests.get(url, **kwargs)


def parse_plddt_from_cif(cif_text: str) -> list[float]:

    plddt_list = []

    cif_lines = cif_text.split("\n")

    cur_residue = -1

    for l in cif_lines:
        cols = l.split()

        if len(cols) == 0 or cols[0] != "ATOM":
            continue

        residue_number = int(cols[8])
        plddt = float(cols[14])

        # Two small sanity checks
        if residue_number < 1:
            raise Exception("Parsed a residue number which is smaller than 1")

        if plddt < 0.0 or plddt > 100.0:
            raise Exception("Parsed a plddt which is not between 1 and 100")

        if residue_number != cur_residue:
            plddt_list.append(plddt)
            cur_residue = residue_number

    return plddt_list


AlphaFoldQueryResult = namedtuple(
    "AlphaFoldQueryResult",
    ["http_status", "accession", "sequence", "plddts", "alpha_fold_data", "cif_text"],
)


def query_alphafold(
    accession: str,
    timeout: int = 10,
    retries: int = 2,
    get_cif: bool = True,
    backoff_time: int = 5,
) -> AlphaFoldQueryResult:
    """For a single accession, query the alpha fold database and retrieve some information."""

    url = ALPHAFOLD_PREDICTION_URL.format(accession=accession)

    # These are the return values
    response_code_alpha_fold: int = -1
    sequence: str | None = None
    plddts: list[float] | None = None
    alpha_fold_data: dict[str, Any] | None = None
    cif_text: str | None = None

    r = get_with_retries(url, retries=retries, backoff_time=5, timeout=timeout)

    response_code_alpha_fold = r.status_code

    if response_code_alpha_fold != 200:
        logger.warning(f"Received response code {response_code_alpha_fold}")
    else:
        try:
            alpha_fold_data = r.json()[0]
        except Exception as e:
            logger.exception("Could not convert response to a single json.")

        if alpha_fold_data is None:
            logger.warning("Received `None` for alpha_fold_data")
        else:
            sequence = alpha_fold_data.get("sequence")

            if sequence is None:
                logger.warning(
                    f"Could not get `sequence` key from alpha fold data. Available keys: {alpha_fold_data.keys()}"
                )
            elif get_cif:
                cif_url = alpha_fold_data.get("cifUrl")

                if cif_url is None:
                    logger.warning(
                        f"Could not get `cifUrl` key from alpha fold data. Available keys: {alpha_fold_data.keys()}"
                    )
                else:

                    r_cif = get_with_retries(
                        cif_url,
                        retries=retries,
                        backoff_time=backoff_time,
                        timeout=timeout,
                    )

                    r_cif.status_code
                    if r_cif.status_code != 200:
                        logger.warning(
                            f"Received response code {r_cif.status_code} for request of cif file"
                        )
                    else:
                        cif_text = r_cif.text
                        plddts = parse_plddt_from_cif(cif_text)
                        if len(sequence) != len(plddts):
                            raise Exception(
                                "sequence and plddts do not have the same length"
                            )

    return AlphaFoldQueryResult(
        response_code_alpha_fold, accession, sequence, plddts, alpha_fold_data, cif_text
    )


def query_alphafold_bulk(
    accession_list: Sequence[str], **kwargs
) -> Generator[AlphaFoldQueryResult]:
    """For a sequence of accessions, return a generator to query the alpha fold database."""

    for a in accession_list:
        yield query_alphafold(a, **kwargs)


if __name__ == "__main__":
    id = "A0A096LP49"
    res = query_alphafold(id)

    print(res.sequence)
    print(res.plddts)
    print(res.cif_text)
    print(res.alpha_fold_data)
    print(res.http_status)

    id_list = [
        "A0A024RBG1",
        "A0A075B6T7",
        "A0A087WTH1",
        "A0A087WTH5",
    ]

    res = query_alphafold_bulk(id_list, get_cif=False)
    print([r.sequence for r in res])
