
germline_snv_probabilities = {
    'A':0.25,
    'C':0.25,
    'T':0.25,
    'G':0.25
}
germline_indel_probabilities = {
    'insertion':0.5,
    'deletion':0.5,
}
somatic_snv_probabilities = {
    'A':0.25,
    'C':0.25,
    'T':0.25,
    'G':0.25
}
somatic_indel_probabilities = {
    'insertion':0.5,
    'deletion':0.5,
}
structural_variations_probabilities = {
    'deletion': 0.2,
    'translocation': 0.2,
    'duplication' : 0.2,
    'inversion' : 0.2,
    'insertion' : 0.2
}


if sum(germline_snv_probabilities.values()) !=1:
    raise ValueError("Error in probabilities_config.py: sum of germline SNV probabilites must equal 1")

if sum(germline_indel_probabilities.values()) !=1:
    raise ValueError("Error in probabilities_config.py: sum of germline InDel probabilites must equal 1")

if sum(germline_snv_probabilities.values()) !=1:
    raise ValueError("Error in probabilities_config.py: sum of somatic SNV probabilites must equal 1")

if sum(germline_snv_probabilities.values()) !=1:
    raise ValueError("Error in probabilities_config.py: sum of somatic InDel probabilites must equal 1")

if sum(germline_snv_probabilities.values()) !=1:
    raise ValueError("Error in probabilities_config.py: sum of somatic tumor SV probabilites must equal 1")
