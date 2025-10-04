from typing import Generator, Any
import requests
import time
from collections.abc import Sequence
from collections import namedtuple

# AlphaFold DB API: returns JSON metadata with model file URLs keyed by UniProt accession
# e.g. https://alphafold.ebi.ac.uk/api/prediction/P69905
ALPHAFOLD_PREDICTION_URL = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"


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
        assert residue_number >= 1
        assert plddt >= 0.0 and plddt <= 100.0

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
    """For a single accession, query the alpha fold database and retrieve the sequence, the plddts, the meta_data and optionally the cif text."""

    url = ALPHAFOLD_PREDICTION_URL.format(accession=accession)

    # These are the return values
    response_code_alpha_fold: int = -1
    sequence: str | None = None
    plddts: list[float] | None = None
    alpha_fold_data: dict[str, Any] | None = None
    cif_text: str | None = None

    response_code_alpha_fold_cif: int | None = None

    for _ in range(retries):
        r = requests.get(url, timeout=timeout)
        response_code_alpha_fold = r.status_code

        if response_code_alpha_fold == 200:
            alpha_fold_data = r.json()[0]
            sequence = alpha_fold_data["sequence"]

            cif_url = alpha_fold_data["cifUrl"]

            if get_cif:
                r_cif = requests.get(cif_url, timeout=timeout)
                response_code_alpha_fold_cif = r_cif.status_code

                if response_code_alpha_fold_cif == 200:
                    cif_text = r_cif.text

                    plddts = parse_plddt_from_cif(cif_text)

                    assert len(sequence) == len(plddts)

                    break
            else:
                break

        time.sleep(backoff_time)

    return AlphaFoldQueryResult(
        response_code_alpha_fold, accession, sequence, plddts, alpha_fold_data, cif_text
    )


def query_alphafold_bulk(
    accession_list: Sequence[str], **kwargs
) -> Generator[AlphaFoldQueryResult]:
    """For a sequence of accession, return a generator query the alpha fold database and retrieve the sequence and the plddts."""

    for a in accession_list:
        yield query_alphafold(a, **kwargs)


if __name__ == "__main__":
    id = "A0A096LP49"
    res = query_alphafold(id)

    print(res.sequence)
    print(res.plddts)
    print(res.cif_text)

    # print(res.alpha_fold_data)
    # print(res.http_status)

    # print(len(seq))
    # print(len(plddts))

    id_list = [
        "A0A024RBG1",
        "A0A075B6T7",
        "A0A087WTH1",
        "A0A087WTH5",
    ]

    res = list(query_alphafold_bulk(id_list))
    print([r.sequence for r in res])
