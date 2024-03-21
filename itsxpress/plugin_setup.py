from q2_types.per_sample_sequences import (SequencesWithQuality,
                                           PairedEndSequencesWithQuality,
                                           JoinedSequencesWithQuality)
from q2_types.sample_data import SampleData
from qiime2.plugin import (Plugin,
                           Str,
                           Choices,
                           Int,
                           Float,
                           Range,
                           Bool,
                           Citations)

from itsxpress.q2_itsxpress import (trim_single,
                                     trim_pair,
                                     trim_pair_output_unmerged,
                                     default_cluster_id)

plugin = Plugin(
    name='itsxpress',
    version='2.0.2',
    package='itsxpress',
    website='https://github.com/USDA-ARS-GBRU/q2_itsxpress             '
            'ITSxpress: https://github.com/USDA-ARS-GBRU/itsxpress',
    description='ITSxpress trims amplicon reads down to their ITS region. '
                'ITSxpress is designed to support the calling of exact sequence variants '
                'rather than OTUs. This newer method of sequence error-correction requires '
                'quality score data from each sequence, so each input sequence must be trimmed. '
                'ITSxpress makes this possible by taking FASTQ data, de-replicating the '
                'sequences then identifying the start and stop sites using HMMSearch. '
                'Results are parsed and the trimmed files are returned. '
                'The ITS 1, ITS2 or the entire ITS region including the 5.8s rRNA'
                'gene can be selected. ALL requires very long reads so it is not routinely'
                'used with Illumina data. ITSxpress uses the hmm models from ITSx so results are comparable.',
    short_description='Plugin for using ITSxpress to rapidly trim the\n'
                      'internally transcribed spacer (ITS) region of FASTQ files.',
    citations=Citations.load('citations.bib', package='itsxpress')
)

taxaList = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'L', 'M', 'ALL', 'O', 'P', 'Q', 'R', 'S', 'T', 'U']

plugin.methods.register_function(
    function=trim_single,
    inputs={'per_sample_sequences': SampleData[SequencesWithQuality]},
    parameters={'region': Str % Choices(['ITS2', 'ITS1', 'ALL']),
                'taxa': Str % Choices(taxaList),
                'threads': Int,
                'cluster_id': Float % Range(0.995, 1.0, inclusive_start=True, inclusive_end=True)},
    outputs=[('trimmed', SampleData[SequencesWithQuality])],
    input_descriptions={'per_sample_sequences': 'The artifact that contains the sequence file(s).'
                                                ' Either Joined Paired or just a single fastq.'
                                                ' One file sequences in the qza data folder.'},
    parameter_descriptions={
        'region': ('\nThe regions ITS2, ITS1, and ALL that can be selected from.'),
        'taxa': ('\nThe selected taxonomic group sequenced that can be selected from.'),
        'threads': ('\nThe number of processor threads to use in the run.'),\
        'cluster_id': ('\nThe percent identity for clustering reads, set to 1 for exact dereplication.')
    },
    output_descriptions={'trimmed': 'The trimmed sequences from ITSxpress.'},
    name='Trim single-end reads',
    description='ITSxpress trimSingle is used for qza types with\n'
                'SquencesWithQuality or JoinedSequencesWithQuality.'
                ' This means the qza must be in the\n'
                'SingleLanePerSampleSingleEndFastqDirFmt, and only contain\n'
                'one fastq file.\n'
                '\nThe taxa list and variable for it:\n'
                '\nA = Alveolata\n'
                '\nB = Bryophyta\n'
                '\nC = Bacillariophyta\n'
                '\nD = Amoebozoa\n'
                '\nE = Euglenozoa\n'
                '\nF = Fungi\n'
                '\nG = Chlorophyta (green algae)\n'
                '\nH = Rhodophyta (red algae)\n'
                '\nI = Phaeophyceae (brown algae)\n'
                '\nL = Marchantiophyta (liverworts)\n'
                '\nM = Metazoa\n'
                '\nO = Oomycota\n'
                '\nP = Haptophyceae (prymnesiophytes)\n'
                '\nQ = Raphidophyceae\n'
                '\nR = Rhizaria\n'
                '\nS = Synurophyceae\n'
                '\nT = Tracheophyta (higher plants)\n'
                '\nU = Eustigmatophyceae\n'
                '\nALL = All'
)

plugin.methods.register_function(
    function=trim_pair,
    inputs={'per_sample_sequences': SampleData[PairedEndSequencesWithQuality]},
    parameters={'region': Str % Choices(['ITS2', 'ITS1', 'ALL']),
                'taxa': Str % Choices(taxaList),
                'threads': Int,
                'reversed_primers': Bool,
                'allow_staggered_reads': Bool,
                'cluster_id': Float % Range(0.995, 1.0, inclusive_start=True, inclusive_end=True)},
    outputs=[('trimmed', SampleData[JoinedSequencesWithQuality])],
    input_descriptions={'per_sample_sequences': 'The artifact that contains the sequence file(s). '
                                                'Only Paired can be used. '
                                                'Two files sequences in the qza data folder'},
    parameter_descriptions={
        'region': ('\nThe regions ITS2, ITS1, and ALL that can be selected from.'),
        'taxa': ('\nThe selected taxonomic group sequenced that can be selected from.'),
        'threads': ('\nThe number of processor threads to use in the run.'),
        'cluster_id': ('\nThe percent identity for clustering reads, set to 1 for exact dereplication.'),
        'allow_staggered_reads': ('\nAllows merging of staggered reads.'),
        'reversed_primers': ('\n Primers are in reverse orientation as in Taylor et al. 2016, DOI:10.1128/AEM.02576-16.')
    },
    output_descriptions={'trimmed': 'The resulting trimmed sequences from ITSxpress'},
    name='Trim paired-end reads, output merged reads for use with Deblur',
    description='ITSxpress trimPair takes the qza type \n'
                'PairedEndSquencesWithQuality. The qza\n'
                'format must be the SingleLanePerSamplePairedEndFastqDirFmt\n'
                'and only contain two fastq files.\n'
                'The function returns merged, trimmed reads for use by Deblur\n'
                '\nThe taxa list and variable for it:\n'
                '\nA = Alveolata\n'
                '\nB = Bryophyta\n'
                '\nC = Bacillariophyta\n'
                '\nD = Amoebozoa\n'
                '\nE = Euglenozoa\n'
                '\nF = Fungi\n'
                '\nG = Chlorophyta (green algae)\n'
                '\nH = Rhodophyta (red algae)\n'
                '\nI = Phaeophyceae (brown algae)\n'
                '\nL = Marchantiophyta (liverworts)\n'
                '\nM = Metazoa\n'
                '\nO = Oomycota\n'
                '\nP = Haptophyceae (prymnesiophytes)\n'
                '\nQ = Raphidophyceae\n'
                '\nR = Rhizaria\n'
                '\nS = Synurophyceae\n'
                '\nT = Tracheophyta (higher plants)\n'
                '\nU = Eustigmatophyceae\n'
                '\nALL = All'

)

plugin.methods.register_function(
    function=trim_pair_output_unmerged,
    inputs={'per_sample_sequences': SampleData[PairedEndSequencesWithQuality]},
    parameters={'region': Str % Choices(['ITS2', 'ITS1', 'ALL']),
                'taxa': Str % Choices(taxaList),
                'threads': Int,
                'reversed_primers': Bool,
                'allow_staggered_reads': Bool,
                'cluster_id': Float % Range(0.995, 1.0, inclusive_start=True, inclusive_end=True)},
    outputs=[('trimmed', SampleData[PairedEndSequencesWithQuality])],
    input_descriptions={'per_sample_sequences': 'The artifact that contains the sequence file(s). '
                                                'Only Paired can be used. '
                                                'Two files sequences in the qza data folder'},
    parameter_descriptions={
        'region': ('\nThe regions ITS2, ITS1, and ALL that can be selected from.'),
        'taxa': ('\nThe selected taxonomic group sequenced that can be selected from.'),
        'threads': ('\nThe number of processor threads to use in the run.'),
        'cluster_id': ('\nThe percent identity for clustering reads, set to 1 for exact dereplication.'),
        'allow_staggered_reads': ('\nAllows merging of staggered reads.'),
        'reversed_primers': ('\n Primers are in reverse orientation as in Taylor et al. 2016, DOI:10.1128/AEM.02576-16.')
    },
    output_descriptions={'trimmed': 'The resulting trimmed sequences from ITSxpress'},
    name='Trim paired-end reads, output unmerged reads for use with Dada2',
    description='ITSxpress trimPairUnmerged takes the qza type \n'
                'PairedEndSquencesWithQuality. The qza\n'
                'format must be the SingleLanePerSamplePairedEndFastqDirFmt\n'
                'and only contain two fastq files.\n'
                'The function returns unmerged, trimmed reads for use by Dada2.\n'
                '\nThe taxa list and variable for it:\n'
                '\nA = Alveolata\n'
                '\nB = Bryophyta\n'
                '\nC = Bacillariophyta\n'
                '\nD = Amoebozoa\n'
                '\nE = Euglenozoa\n'
                '\nF = Fungi\n'
                '\nG = Chlorophyta (green algae)\n'
                '\nH = Rhodophyta (red algae)\n'
                '\nI = Phaeophyceae (brown algae)\n'
                '\nL = Marchantiophyta (liverworts)\n'
                '\nM = Metazoa\n'
                '\nO = Oomycota\n'
                '\nP = Haptophyceae (prymnesiophytes)\n'
                '\nQ = Raphidophyceae\n'
                '\nR = Rhizaria\n'
                '\nS = Synurophyceae\n'
                '\nT = Tracheophyta (higher plants)\n'
                '\nU = Eustigmatophyceae\n'
                '\nALL = All'

)
