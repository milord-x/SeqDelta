from seqdelta.mutation import analyze_sequences

REFERENCE = "ATGGAATTTCCGTGGAAATAA"


def test_substitution_detection():
    result = analyze_sequences(REFERENCE, "ATGGACTTTCCGTGGAAATAA")
    assert result.summary.total_mutations == 1
    assert result.summary.substitutions == 1
    mutation = result.nucleotide_mutations[0]
    assert mutation.mutation_type == "substitution"
    assert mutation.position == 6
    assert mutation.ref_nt == "A"
    assert mutation.mut_nt == "C"


def test_insertion_detection():
    result = analyze_sequences(REFERENCE, "ATGGAATTTACCGTGGAAATAA")
    assert result.summary.insertions == 1
    mutation = result.nucleotide_mutations[0]
    assert mutation.mutation_type == "insertion"
    assert mutation.mut_nt == "A"
    assert mutation.effect == "frameshift"


def test_deletion_detection():
    result = analyze_sequences(REFERENCE, "ATGGAATTTCGTGGAAATAA")
    assert result.summary.deletions == 1
    mutation = result.nucleotide_mutations[0]
    assert mutation.mutation_type == "deletion"
    assert mutation.ref_nt == "C"
    assert mutation.effect == "frameshift"


def test_silent_mutation():
    result = analyze_sequences(REFERENCE, "ATGGAGTTTCCGTGGAAATAA")
    mutation = result.nucleotide_mutations[0]
    assert mutation.effect == "silent"
    assert mutation.ref_codon == "GAA"
    assert mutation.mut_codon == "GAG"
    assert mutation.ref_aa == "E"
    assert mutation.mut_aa == "E"


def test_missense_mutation():
    result = analyze_sequences(REFERENCE, "ATGGACTTTCCGTGGAAATAA")
    mutation = result.nucleotide_mutations[0]
    assert mutation.effect == "missense"
    assert mutation.ref_aa == "E"
    assert mutation.mut_aa == "D"


def test_nonsense_mutation():
    result = analyze_sequences(REFERENCE, "ATGTAATTTCCGTGGAAATAA")
    mutation = result.nucleotide_mutations[0]
    assert mutation.effect == "nonsense"
    assert mutation.ref_codon == "GAA"
    assert mutation.mut_codon == "TAA"
    assert mutation.mut_aa == "*"


def test_frameshift_summary():
    result = analyze_sequences(REFERENCE, "ATGGAATTTACCGTGGAAATAA")
    assert result.summary.frameshift == 1
    changed_codons = [codon for codon in result.codon_changes if codon.changed]
    assert changed_codons
    assert changed_codons[0].effect in {"missense", "frameshift"}
    downstream = [codon for codon in changed_codons if codon.effect == "frameshift"]
    assert downstream


def test_identity_counts_mismatches():
    result = analyze_sequences(REFERENCE, "ATGGACTTTCCGTGGAAATAA")
    assert result.summary.identity < 1.0
    assert round(result.summary.identity * 100, 1) == 95.2
