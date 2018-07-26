
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

### TO DO! (but arguably not necessary)
### we should probably allow users to set average tumor SV lengths via
### probabilities_config.py
### prob = np.random.uniform(0.001, 0.0000001, 1)   ## draw prob from uniform, 0.001 to 1e-7
### some silly parameter to "tune" if they wanted

### same with duplications parameters
### duplication_prob = np.random.uniform(0.05, 0.7, 1)  ## with np.random.geometric(p, 1), these values are 14 to 2
## num_duplications = self.get_event_length(p=duplication_prob[0]) # exponential ranging from 1 to 10

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
