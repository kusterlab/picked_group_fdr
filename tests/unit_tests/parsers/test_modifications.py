import picked_group_fdr.parsers.modifications as mod


class TestMaxQuantToProForma:
    def test_maxquant_to_internal_carbamidomethylation(self):
        assert mod.maxquant_mod_to_proforma()("_ABCDEFGH_") == "ABC[UNIMOD:4]DEFGH"

    def test_maxquant_to_internal_tmt(self):
        fixed_mods = {"C": "C[UNIMOD:4]", "^_": "_[UNIMOD:737]-", "K": "K[UNIMOD:737]"}
        assert (
            mod.maxquant_mod_to_proforma(fixed_mods)("_ABCDEFGHK_")
            == "[UNIMOD:737]-ABC[UNIMOD:4]DEFGHK[UNIMOD:737]"
        )

    def test_maxquant_to_internal_silac(self):
        fixed_mods = {"C": "C[UNIMOD:4]", "K": "K[UNIMOD:259]", "R": "R[UNIMOD:267]"}
        assert (
            mod.maxquant_mod_to_proforma(fixed_mods)("_ABCDEFGHRK_")
            == "ABC[UNIMOD:4]DEFGHR[UNIMOD:267]K[UNIMOD:259]"
        )


class TestPrositToProForma:
    def test_prosit_to_internal_carbamidomethylation(self):
        assert (
            mod.prosit_mod_to_proforma()("ABC[UNIMOD:4]DEFGH") == "ABC[UNIMOD:4]DEFGH"
        )

    def test_prosit_to_internal_tmt(self):
        assert (
            mod.prosit_mod_to_proforma()("[UNIMOD:737]ABC[UNIMOD:4]DEFGHK[UNIMOD:737]")
            == "[UNIMOD:737]-ABC[UNIMOD:4]DEFGHK[UNIMOD:737]"
        )

    def test_prosit_c_terminal_mod(self):
        assert (
            mod.prosit_mod_to_proforma(
                variable_mods={"[UNIMOD:737]$": "-[UNIMOD:737]"}
            )("[UNIMOD:737]ABC[UNIMOD:4]DEFGHK[UNIMOD:737]")
            == "[UNIMOD:737]-ABC[UNIMOD:4]DEFGHK-[UNIMOD:737]"
        )
