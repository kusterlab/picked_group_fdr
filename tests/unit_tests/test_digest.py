from unittest.mock import mock_open, patch

import picked_group_fdr.digest as digest
from picked_group_fdr.digestion_params import DigestionParams


class TestIsValidPrositPeptide:
    def test_valid_peptide(self):
        assert digest.is_valid_prosit_peptide("ACDEFGHIKLMNPQRSTVWY")
        assert digest.is_valid_prosit_peptide("AAAAA")

    def test_invalid_length(self):
        # Peptide length exceeds 30 characters
        assert not digest.is_valid_prosit_peptide("ACDEFGHIKLMNPQRSTVWYACDEFGHIKLM")

    def test_invalid_characters(self):
        # Peptide contains 'U'
        assert not digest.is_valid_prosit_peptide("ACDEFGHIUKLMNPQRSTVWY")

        # Peptide contains 'X'
        assert not digest.is_valid_prosit_peptide("AXDEFGHIKLMNPQRSTVWY")

    def test_invalid_length_and_characters(self):
        # Peptide contains 'U' and exceeds 30 characters
        assert not digest.is_valid_prosit_peptide("ACDEFGHIKLMNPQRSTVWYU" * 2)

        # Peptide contains 'X' and exceeds 30 characters
        assert not digest.is_valid_prosit_peptide("AXDEFGHIKLMNPQRSTVWY" * 2)


class TestParseUntilFirstSpace:
    def test_valid_fasta_id(self):
        fasta_id = "protein123 description"
        result = digest.parse_until_first_space(fasta_id)
        assert result == "protein123"

    def test_fasta_id_without_space(self):
        fasta_id = "protein123"
        result = digest.parse_until_first_space(fasta_id)
        assert result == "protein123"

    def test_empty_fasta_id(self):
        fasta_id = ""
        result = digest.parse_until_first_space(fasta_id)
        assert result == ""

    def test_fasta_id_with_multiple_spaces(self):
        fasta_id = "protein123 description with spaces"
        result = digest.parse_until_first_space(fasta_id)
        assert result == "protein123"

    def test_fasta_id_with_leading_space(self):
        fasta_id = " protein123 description"
        result = digest.parse_until_first_space(fasta_id)
        assert result == ""

    def test_fasta_id_with_trailing_space(self):
        fasta_id = "protein123 description "
        result = digest.parse_until_first_space(fasta_id)
        assert result == "protein123"


def test_swap_special_aas():
    seq = "ABCKDEFRRR"
    specialAAs = ["K", "R"]
    assert digest.swap_special_aas(seq, specialAAs) == "ABKCDERRRF"


def test_has_miscleavage_true():
    seq = "ABCDEFKR"
    assert digest.has_miscleavage(seq) == True


def test_has_miscleavage_false():
    seq = "ABCDEFKPR"
    assert digest.has_miscleavage(seq) == False


class TestNonSpecificDigest:
    def test_non_specific_digest(self):
        seq = "ABCDEFGH"
        min_len = 6
        max_len = 30
        assert set(digest.non_specific_digest(seq, min_len, max_len)) == set(
            ["ABCDEF", "BCDEFG", "CDEFGH", "ABCDEFG", "BCDEFGH", "ABCDEFGH"]
        )

    def test_non_specific_digest_max_len(self):
        seq = "ABCDEFGH"
        min_len = 6
        max_len = 7
        assert set(digest.non_specific_digest(seq, min_len, max_len)) == set(
            ["ABCDEF", "BCDEFG", "CDEFGH", "ABCDEFG", "BCDEFGH"]
        )


class TestSemiSpecificDigest:
    def test_semi_specific_digest_no_cleavage_site(self):
        seq = "ABCDEFGH"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 2
        methionineCleavage = True
        assert set(
            digest.semi_specific_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(["ABCDEF", "ABCDEFG", "ABCDEFGH", "BCDEFGH", "CDEFGH"])

    def test_semi_specific_digest_methionine_cleavage(self):
        seq = "MABCDEFGH"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 2
        methionineCleavage = True
        assert set(
            digest.semi_specific_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(
            [
                "MABCDE",
                "MABCDEF",
                "MABCDEFG",
                "MABCDEFGH",
                "ABCDEF",
                "ABCDEFG",
                "ABCDEFGH",
                "BCDEFGH",
                "CDEFGH",
            ]
        )

    # make sure that the methionine cleavage is not counted as a miscleavage
    def test_semi_specific_digest_methionine_cleavage_plus_one_miscleavage(self):
        seq = "MABCDEFKKK"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 1
        methionineCleavage = True
        assert set(
            digest.semi_specific_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(
            [
                "MABCDE",
                "MABCDEF",
                "MABCDEFK",
                "MABCDEFKK",
                "ABCDEF",
                "ABCDEFK",
                "ABCDEFKK",
                "BCDEFK",
                "BCDEFKK",
                "CDEFKK",
            ]
        )

    def test_semi_specific_digest_no_cleavage_site_max_len(self):
        seq = "ABCDEFGH"
        min_len = 6
        max_len = 7
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 2
        methionineCleavage = True
        assert set(
            digest.semi_specific_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(["ABCDEF", "ABCDEFG", "BCDEFGH", "CDEFGH"])

    def test_semi_specific_digest_no_miscleavage(self):
        seq = "ABCDEFGKX"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 0
        methionineCleavage = True
        assert set(
            digest.semi_specific_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(["ABCDEF", "ABCDEFG", "ABCDEFGK", "BCDEFGK", "CDEFGK"])

    def test_semi_specific_digest_one_miscleavage(self):
        seq = "ABCDEFGKX"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 1
        methionineCleavage = True
        assert set(
            digest.semi_specific_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(
            [
                "ABCDEF",
                "ABCDEFG",
                "ABCDEFGK",
                "BCDEFGK",
                "CDEFGK",
                "ABCDEFGKX",
                "BCDEFGKX",
                "CDEFGKX",
                "DEFGKX",
            ]
        )

    def test_semi_specific_digest_not_post(self):
        seq = "ABCDEFKPX"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 0
        methionineCleavage = True
        assert set(
            digest.semi_specific_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(
            [
                "ABCDEF",
                "ABCDEFK",
                "ABCDEFKP",
                "ABCDEFKPX",
                "BCDEFKPX",
                "CDEFKPX",
                "DEFKPX",
            ]
        )


class TestFullDigest:
    def test_full_digest_no_cleavage_site(self):
        seq = "ABCDEFGH"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 2
        methionineCleavage = True
        assert set(
            digest.full_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(["ABCDEFGH"])

    def test_full_digest_methionine_cleavage(self):
        seq = "MABCDEFGH"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 2
        methionineCleavage = True
        assert set(
            digest.full_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(["MABCDEFGH", "ABCDEFGH"])

    # make sure that the methionine cleavage is not counted as a miscleavage
    def test_full_digest_methionine_cleavage_one_miscleavage(self):
        seq = "MABCDEFGHKKK"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 1
        methionineCleavage = True
        assert set(
            digest.full_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(["MABCDEFGHK", "MABCDEFGHKK", "ABCDEFGHK", "ABCDEFGHKK"])

    def test_full_digest_no_cleavage_site_max_len(self):
        seq = "ABCDEFGH"
        min_len = 6
        max_len = 7
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 2
        methionineCleavage = True
        assert set(
            digest.full_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set([])

    def test_full_digest_no_miscleavage(self):
        seq = "ABCDEFGKX"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 0
        methionineCleavage = True
        assert set(
            digest.full_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(["ABCDEFGK"])

    def test_full_digest_one_miscleavage(self):
        seq = "ABCDEFGKX"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 1
        methionineCleavage = True
        assert set(
            digest.full_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(["ABCDEFGK", "ABCDEFGKX"])

    def test_full_digest_not_post(self):
        seq = "ABCDEFKPX"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 0
        methionineCleavage = True
        assert set(
            digest.full_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(["ABCDEFKPX"])

    def test_full_digest_post(self):
        seq = "ABCDEFKPXAAA"
        min_len = 6
        max_len = 30
        pre = [""]
        not_post = [""]
        post = ["K"]
        miscleavages = 0
        methionineCleavage = True
        assert set(
            digest.full_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(["ABCDEF", "KPXAAA"])


class TestReadFastaMaxQuant:
    @staticmethod
    def _run_read_fasta_maxquant(
        file_content,
        db="target",
        parse_id=digest.parse_until_first_space,
        special_aas=None,
    ):
        with patch("builtins.open", mock_open(read_data=file_content)):
            return list(
                digest.read_fasta_maxquant(
                    "test.fasta", db=db, parse_id=parse_id, special_aas=special_aas
                )
            )

    def test_target_db(self):
        fasta_content = ">protein1\nSEQUENCE1\n>protein2\nSEQUENCE2\n"
        result = self._run_read_fasta_maxquant(fasta_content, db="target")
        assert result == [("protein1", "SEQUENCE1"), ("protein2", "SEQUENCE2")]

    def test_decoy_db(self):
        fasta_content = ">protein1\nSEQUENCE1\n>protein2\nSEQUENCE2\n"
        result = self._run_read_fasta_maxquant(fasta_content, db="decoy")
        assert result == [
            ("REV__protein1", "1ECNEUQES"),
            ("REV__protein2", "2ECNEUQES"),
        ]

    def test_concat_db(self):
        fasta_content = ">protein1\nSEQUENCE1\n>protein2\nSEQUENCE2\n"
        result = self._run_read_fasta_maxquant(fasta_content, db="concat")
        assert result == [
            ("protein1", "SEQUENCE1"),
            ("REV__protein1", "1ECNEUQES"),
            ("protein2", "SEQUENCE2"),
            ("REV__protein2", "2ECNEUQES"),
        ]

    def test_custom_parse_id(self):
        def custom_parse_id(identifier: str) -> str:
            return identifier.upper()

        fasta_content = ">protein1\nSEQUENCE1\n>protein2\nSEQUENCE2\n"
        result = self._run_read_fasta_maxquant(fasta_content, parse_id=custom_parse_id)
        assert result == [("PROTEIN1", "SEQUENCE1"), ("PROTEIN2", "SEQUENCE2")]

    def test_with_special_aas(self):
        fasta_content = ">protein1\nSEQUENCE1\n>protein2\nSEQUENCE2\n"
        result = self._run_read_fasta_maxquant(
            fasta_content, db="concat", special_aas=["N"]
        )
        assert result == [
            ("protein1", "SEQUENCE1"),
            ("REV__protein1", "1ENCEUQES"),
            ("protein2", "SEQUENCE2"),
            ("REV__protein2", "2ENCEUQES"),
        ]


class TestGetProteinSequences:
    @staticmethod
    def _run_get_protein_sequences(file_contents, **kwargs):
        with patch("builtins.open", mock_open(read_data=file_contents)):
            return digest.get_protein_sequences(["test.fasta"], **kwargs)

    def test_empty_file_paths(self):
        result = digest.get_protein_sequences(None)
        assert result == {}

    def test_empty_file(self):
        result = self._run_get_protein_sequences(None)
        assert result == {}

    def test_single_file(self):
        fasta_content = ">protein1\nSEQUENCE1\n>protein2\nSEQUENCE2\n"
        result = self._run_get_protein_sequences(fasta_content)
        assert result == {"protein1": "SEQUENCE1", "protein2": "SEQUENCE2"}

    def test_multiple_files(self):
        fasta_content1 = ">protein1\nSEQUENCE1\n>protein2\nSEQUENCE2\n"
        fasta_content2 = ">protein3\nSEQUENCE3\n>protein4\nSEQUENCE4\n"
        result = self._run_get_protein_sequences(fasta_content1 + fasta_content2)
        assert result == {
            "protein1": "SEQUENCE1",
            "protein2": "SEQUENCE2",
            "protein3": "SEQUENCE3",
            "protein4": "SEQUENCE4",
        }

    def test_duplicate_protein_ids(self):
        fasta_content = ">protein1\nSEQUENCE1\n>protein1\nSEQUENCE2\n"
        result = self._run_get_protein_sequences(fasta_content)
        assert result == {"protein1": "SEQUENCE1"}

    def test_custom_read_fasta_kwargs(self):
        fasta_content = ">protein1\nSEQUENCE1\n>protein2\nSEQUENCE2\n"
        result = self._run_get_protein_sequences(
            fasta_content, db="target", parse_id=lambda x: x.upper()
        )
        assert result == {"PROTEIN1": "SEQUENCE1", "PROTEIN2": "SEQUENCE2"}


class TestPeptideProteinMappingFunctions:
    @staticmethod
    def _run_get_peptide_to_protein_map_from_params(
        fasta_content, digestion_params_list
    ):
        with patch("builtins.open", mock_open(read_data=fasta_content)):
            return digest.get_peptide_to_protein_map_from_params(
                ["test.fasta"], digestion_params_list
            )

    @staticmethod
    def _run_get_peptide_to_protein_map_from_params_single(
        fasta_content, digestion_params
    ):
        with patch("builtins.open", mock_open(read_data=fasta_content)):
            return digest.get_peptide_to_protein_map_from_params_single(
                "test.fasta", digestion_params
            )

    @staticmethod
    def _run_get_peptide_to_protein_map(fasta_content, **kwargs):
        with patch("builtins.open", mock_open(read_data=fasta_content)):
            return digest.get_peptide_to_protein_map("test.fasta", **kwargs)

    def test_get_peptide_to_protein_map_from_params(self):
        fasta_content = ">protein1\nSEQUENCE1\n>protein2\nSEQUENCE2\n"
        digestion_params_list = [DigestionParams(enzyme="trypsin", digestion="full")]
        result = self._run_get_peptide_to_protein_map_from_params(
            fasta_content, digestion_params_list
        )
        assert result == {
            "SEQUENCE1": ["protein1"],
            "SEQUENCE2": ["protein2"],
            "1ECNEUQES": ["REV__protein1"],
            "2ECNEUQES": ["REV__protein2"],
        }

    def test_get_peptide_to_protein_map_from_params_single(self):
        fasta_content = ">protein1\nSEQUENCE1\n>protein2\nSEQUENCE2\n"
        digestion_params = DigestionParams(enzyme="trypsin", digestion="full")
        result = self._run_get_peptide_to_protein_map_from_params_single(
            fasta_content, digestion_params
        )
        assert result == {
            "SEQUENCE1": ["protein1"],
            "SEQUENCE2": ["protein2"],
            "1ECNEUQES": ["REV__protein1"],
            "2ECNEUQES": ["REV__protein2"],
        }

    def test_get_peptide_to_protein_map(self):
        fasta_content = ">protein1\nSEQUENCE1\n>protein2\nSEQUENCE2\n"
        result = self._run_get_peptide_to_protein_map(fasta_content)
        assert result == {
            "SEQUENCE1": ["protein1"],
            "SEQUENCE2": ["protein2"],
            "1ECNEUQES": ["REV__protein1"],
            "2ECNEUQES": ["REV__protein2"],
        }

    def test_get_peptide_to_protein_map_with_use_hash_key_true(self):
        fasta_content = ">protein1\nSEQUENCE1\n>protein2\nSEQUENCE2\n"
        kwargs = {
            "db": "concat",
            "min_len": 5,
            "max_len": 10,
            "pre": ["K"],
            "not_post": ["P"],
            "post": ["S"],
            "digestion": "partial",
            "miscleavages": 1,
            "methionine_cleavage": False,
            "use_hash_key": True,
            "special_aas": ["K"],
            "parse_id": lambda x: x.upper(),
        }
        result = self._run_get_peptide_to_protein_map(fasta_content, **kwargs)
        assert result == (
            {
                "SEQUEN": ["PROTEIN1", "PROTEIN2"],
                "1ECNEU": ["REV__PROTEIN1"],
                "2ECNEU": ["REV__PROTEIN2"],
            },
            {
                "PROTEIN1": "SEQUENCE1",
                "PROTEIN2": "SEQUENCE2",
                "REV__PROTEIN1": "1ECNEUQES",
                "REV__PROTEIN2": "2ECNEUQES",
            },
        )

    def test_get_peptide_to_protein_map_with_custom_args(self):
        fasta_content = ">protein1\nSEQUENCE1\n>protein2\nSEQUENCE2\n"
        kwargs = {
            "db": "concat",
            "min_len": 5,
            "max_len": 10,
            "pre": ["K"],
            "not_post": ["P"],
            "post": ["S"],
            "digestion": "partial",
            "miscleavages": 1,
            "methionine_cleavage": False,
            "use_hash_key": False,
            "special_aas": ["K"],
            "parse_id": lambda x: x.upper(),
        }
        result = self._run_get_peptide_to_protein_map(fasta_content, **kwargs)
        assert result == {
            "SEQUENCE1": ["PROTEIN1"],
            "SEQUENCE2": ["PROTEIN2"],
            "1ECNEUQES": ["REV__PROTEIN1"],
            "2ECNEUQES": ["REV__PROTEIN2"],
            "1ECNEUQE": ["REV__PROTEIN1"],
            "2ECNEUQE": ["REV__PROTEIN2"],
        }

    def test_get_peptide_to_protein_map_with_use_hash_key_false(self):
        fasta_content = ">protein1\nSEQUENCE1\n>protein2\nSEQUENCE2\n"
        result = self._run_get_peptide_to_protein_map(fasta_content, use_hash_key=False)
        assert result == {
            "SEQUENCE1": ["protein1"],
            "SEQUENCE2": ["protein2"],
            "1ECNEUQES": ["REV__protein1"],
            "2ECNEUQES": ["REV__protein2"],
        }
