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
                           Citations,
                           Collection)

from itsxpress.q2_itsxpress import (trim_single,
                                     trim_pair,
                                     trim_pair_output_unmerged,
                                     split_single_end,
                                     split_paired_end,
                                     combine_single,
                                     combine_pair,
                                     combine_pair_unmerged,
                                     default_cluster_id)
from itsxpress.pipelines import (trim_single as trim_single_pipeline,
                                  trim_pair as trim_pair_pipeline,
                                  trim_pair_output_unmerged as trim_pair_output_unmerged_pipeline)
from ._version import __version__

trim_single_sample = trim_single
trim_single_sample.__name__ = "trim_single_sample"

trim_pair_sample = trim_pair
trim_pair_sample.__name__ = "trim_pair_sample"

trim_pair_sample_unmerged = trim_pair_output_unmerged
trim_pair_sample_unmerged.__name__ = "trim_pair_sample_unmerged"
plugin = Plugin(
    name='itsxpress',
    version=__version__,
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

taxaList = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'L', 'M', 'ALL', 'O', 'P', 'Q', 'R', 'S', 'T', 'U','Y']

plugin.methods.register_function(
    function=trim_single_sample,
    inputs={'per_sample_sequences': SampleData[SequencesWithQuality]},
    parameters={'region': Str % Choices(['ITS2', 'ITS1', 'ALL']),
                'taxa': Str % Choices(taxaList),
                'threads': Int,
                'cluster_id': Float % Range(0.995, 1.0, inclusive_start=True, inclusive_end=True),
                'trim_ccs': Bool},
    outputs=[('trimmed', SampleData[SequencesWithQuality])],
    input_descriptions={'per_sample_sequences': 'The artifact that contains the sequence file(s).'
                                                ' Either Joined Paired or just a single fastq.'
                                                ' One file sequences in the qza data folder.'},
    parameter_descriptions={
        'region': ('\nThe regions ITS2, ITS1, and ALL that can be selected from.'),
        'taxa': ('\nThe selected taxonomic group sequenced that can be selected from.'),
        'threads': ('\nThe number of processor threads to use in the run.'),\
        'cluster_id': ('\nThe percent identity for clustering reads, set to 1 for exact dereplication.'),
        'trim_ccs': ('\nWhether to enable trim-ccs mode for PacBio CCS reads.')
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
                '\nY = Parabasalia\n'
                '\nALL = All'
)

plugin.methods.register_function(
    function=trim_pair_sample,
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
                '\nY = Parabasalia\n'
                '\nALL = All'

)

plugin.methods.register_function(
    function=split_single_end,
    inputs={'per_sample_sequences': SampleData[SequencesWithQuality]},
    parameters={},
    outputs=[('splits', Collection[SampleData[SequencesWithQuality]])],
    input_descriptions={'per_sample_sequences': 'The single-end sequence artifact to split.'},
    parameter_descriptions={},
    output_descriptions={'splits': 'The individual single-end sequence splits.'},
    name='Split single-end reads',
    description='Splits single-end sequence artifact into individual sample artifacts.'
)

plugin.methods.register_function(
    function=split_paired_end,
    inputs={'per_sample_sequences': SampleData[PairedEndSequencesWithQuality]},
    parameters={},
    outputs=[('splits', Collection[SampleData[PairedEndSequencesWithQuality]])],
    input_descriptions={'per_sample_sequences': 'The paired-end sequence artifact to split.'},
    parameter_descriptions={},
    output_descriptions={'splits': 'The individual paired-end sequence splits.'},
    name='Split paired-end reads',
    description='Splits paired-end sequence artifact into individual sample artifacts.'
)

plugin.methods.register_function(
    function=combine_single,
    inputs={'results': Collection[SampleData[SequencesWithQuality]]},
    parameters={},
    outputs=[('combined', SampleData[SequencesWithQuality])],
    input_descriptions={'results': 'The trimmed single-end sequence splits to combine.'},
    parameter_descriptions={},
    output_descriptions={'combined': 'The combined trimmed single-end sequence artifact.'},
    name='Combine single-end reads',
    description='Combines trimmed single-end sequence splits back into one unified artifact.'
)

plugin.methods.register_function(
    function=combine_pair,
    inputs={'results': Collection[SampleData[JoinedSequencesWithQuality]]},
    parameters={},
    outputs=[('combined', SampleData[JoinedSequencesWithQuality])],
    input_descriptions={'results': 'The trimmed merged paired-end sequence splits to combine.'},
    parameter_descriptions={},
    output_descriptions={'combined': 'The combined trimmed merged paired-end sequence artifact.'},
    name='Combine merged paired-end reads',
    description='Combines trimmed merged paired-end sequence splits back into one unified artifact.'
)

plugin.methods.register_function(
    function=combine_pair_unmerged,
    inputs={'results': Collection[SampleData[PairedEndSequencesWithQuality]]},
    parameters={},
    outputs=[('combined', SampleData[PairedEndSequencesWithQuality])],
    input_descriptions={'results': 'The trimmed unmerged paired-end sequence splits to combine.'},
    parameter_descriptions={},
    output_descriptions={'combined': 'The combined trimmed unmerged paired-end sequence artifact.'},
    name='Combine unmerged paired-end reads',
    description='Combines trimmed unmerged paired-end sequence splits back into one unified artifact.'
)


plugin.pipelines.register_function(
    function=trim_single_pipeline,
    inputs={'per_sample_sequences': SampleData[SequencesWithQuality]},
    parameters={'region': Str % Choices(['ITS2', 'ITS1', 'ALL']),
                'taxa': Str % Choices(taxaList),
                'threads': Int,
                'cluster_id': Float % Range(0.995, 1.0, inclusive_start=True, inclusive_end=True),
                'trim_ccs': Bool,
                'num_splits': Int},
    outputs=[('trimmed', SampleData[SequencesWithQuality])],
    input_descriptions={'per_sample_sequences': 'The artifact that contains the sequence file(s).'},
    parameter_descriptions={
        'region': '\nThe regions ITS2, ITS1, and ALL that can be selected from.',
        'taxa': '\nThe selected taxonomic group sequenced that can be selected from.',
        'threads': '\nThe number of processor threads to use in the run.',
        'cluster_id': '\nThe percent identity for clustering reads, set to 1 for exact dereplication.',
        'trim_ccs': '\nWhether to enable trim-ccs mode for PacBio CCS reads.',
        'num_splits': '\nThe number of splits to create (unused/deprecated).'
    },
    output_descriptions={'trimmed': 'The trimmed sequences from ITSxpress.'},
    name='Trim single-end reads',
    description='Pipeline for trimming of single-end reads.'
)

plugin.pipelines.register_function(
    function=trim_pair_pipeline,
    inputs={'per_sample_sequences': SampleData[PairedEndSequencesWithQuality]},
    parameters={'region': Str % Choices(['ITS2', 'ITS1', 'ALL']),
                'taxa': Str % Choices(taxaList),
                'threads': Int,
                'reversed_primers': Bool,
                'allow_staggered_reads': Bool,
                'cluster_id': Float % Range(0.995, 1.0, inclusive_start=True, inclusive_end=True),
                'num_splits': Int},
    outputs=[('trimmed', SampleData[JoinedSequencesWithQuality])],
    input_descriptions={'per_sample_sequences': 'The artifact that contains the sequence file(s).'},
    parameter_descriptions={
        'region': '\nThe regions ITS2, ITS1, and ALL that can be selected from.',
        'taxa': '\nThe selected taxonomic group sequenced that can be selected from.',
        'threads': '\nThe number of processor threads to use in the run.',
        'cluster_id': '\nThe percent identity for clustering reads, set to 1 for exact dereplication.',
        'allow_staggered_reads': '\nAllows merging of staggered reads.',
        'reversed_primers': '\n Primers are in reverse orientation.',
        'num_splits': '\nThe number of splits to create (unused/deprecated).'
    },
    output_descriptions={'trimmed': 'The resulting trimmed sequences from ITSxpress.'},
    name='Trim paired-end reads, output merged reads for use with Deblur',
    description='Pipeline for trimming of paired-end reads (merged output).'
)

plugin.pipelines.register_function(
    function=trim_pair_output_unmerged_pipeline,
    inputs={'per_sample_sequences': SampleData[PairedEndSequencesWithQuality]},
    parameters={'region': Str % Choices(['ITS2', 'ITS1', 'ALL']),
                'taxa': Str % Choices(taxaList),
                'threads': Int,
                'reversed_primers': Bool,
                'allow_staggered_reads': Bool,
                'cluster_id': Float % Range(0.995, 1.0, inclusive_start=True, inclusive_end=True),
                'num_splits': Int},
    outputs=[('trimmed', SampleData[PairedEndSequencesWithQuality])],
    input_descriptions={'per_sample_sequences': 'The artifact that contains the sequence file(s).'},
    parameter_descriptions={
        'region': '\nThe regions ITS2, ITS1, and ALL that can be selected from.',
        'taxa': '\nThe selected taxonomic group sequenced that can be selected from.',
        'threads': '\nThe number of processor threads to use in the run.',
        'cluster_id': '\nThe percent identity for clustering reads, set to 1 for exact dereplication.',
        'allow_staggered_reads': '\nAllows merging of staggered reads.',
        'reversed_primers': '\n Primers are in reverse orientation.',
        'num_splits': '\nThe number of splits to create (unused/deprecated).'
    },
    output_descriptions={'trimmed': 'The resulting trimmed sequences from ITSxpress.'},
    name='Trim paired-end reads, output unmerged reads for use with Dada2',
    description='Pipeline for trimming of paired-end reads (unmerged output).'
)

plugin.methods.register_function(
    function=trim_pair_sample_unmerged,
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
                '\nY = Parabasalia\n'
                '\nALL = All'

)
