# PromoterOntology

The goal of the Promotor Ontology project is to annotate promoters with three sets of features:

    Transcription factor binding sites,
    CAGE features
    Other sequence features (CpG islands, bidirectionality....)

Just like with Gene Ontology where for a given group of genes, we can screen for overrepresented features of given group. With Promoter Ontology we will be able to find overrepresented promoter features in a given group of genes.

For every new tissue analysis, this is the order of steps to take in order to get to Ontology stage:

    Calculate CAGE TCs for tissue/sample of interest
    From TCs create promoter regions 500 bp from each side of dTSS
    link these promoters to genes
    Fetch PWMs for TFs of interest
    Scan promoter sequeneces with obtained PWMs
    Normalise hits for location frequency and score
    threshold all hits
    Create matrix where for each promoter summarise which TFs are most likely bound
    Run hypergeometric test for set of genes provided

