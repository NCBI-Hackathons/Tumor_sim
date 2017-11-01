
structural_variations_probabilities = {
        'deletion': 0.2,
        'translocation': 0.2,
        'duplication' : 0.2,
        'inversion' : 0.2,
        'insertion' : 0.2
}
snv_probabilities = {
    'A':.25,
    'C':.25,
    'T':.25,
    'G':.25
}

variant_distribution = 'uniform'



if variant_distribution not in ('uniform', 'normal', 'gaussian'):
    raise NotImplementedError("Only Uniform and Gaussian are implemented! Enter 'uniform' or 'normal' ")
