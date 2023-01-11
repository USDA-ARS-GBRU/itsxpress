import logging

class ItsPosition:
    """Class for ITS positional information derived from hmmserach domtable files.
    Args:
        domtable (str):  the path locating the domtable file from HMMER 3 hmmsearch.
        region (str): The region of the ITS to extract choices: ["ITS1", "ITS2", "ALL"].
    Attributes:
        ddict (dict): A dictionary holding the scores and start and stop
            positions for the selected segment of each sequence.
            Example: {sample:{left:{score:31, pos:15}, right:{score:32, pos:354}, tlen:449}
        leftprefix (str): the left prefix to search for (set by type variable).
        rightprefix (str): the right prefix to search for (set by type variable).
    """
    def _score(self, sequence, stype, score, from_pos, to_pos, tlen):
        """Evaluates scores and positions from the new line of a domtable file and
            updates ddict if necessary.
        Args:
            sequence (str): The name of the sequence.
            stype (str): {'left', 'right'}
            score (int): The bit score from HMMSearch.
            to_pos (int): the ending position of the left sequence.
            from_pos (int): The beginning position of the right sequence.
            tlen (int): the length of the sequence
        """

        if stype in self.ddict[sequence]:
            if score > self.ddict[sequence][stype]["score"]:
                self.ddict[sequence][stype]["score"] = score
                self.ddict[sequence][stype]["to_pos"] = to_pos
                self.ddict[sequence][stype]["from_pos"] = from_pos
        else:
            self.ddict[sequence][stype] = {}
            self.ddict[sequence][stype]["score"] = score
            self.ddict[sequence][stype]["to_pos"] = to_pos
            self.ddict[sequence][stype]["from_pos"] = from_pos
            self.ddict[sequence]["tlen"] = tlen


    def parse(self):
        """Parses dom table from HMMsearch.
        The dom table is parsed and the start and stop position from the top scoring
        hmm math is saved. The start and stop positions of reach sequence are added to the ddict attribute.
        """
        try:
            with open(self.domtable, 'r') as f:
                for num, line in enumerate(f):
                    if not line.startswith("#"):
                        ll = line.split()
                        sequence = ll[0]
                        hmmprofile = ll[3]
                        score = float(ll[13])
                        from_pos = int(ll[19])
                        to_pos = int(ll[20])
                        tlen = int(ll[2])
                        if sequence not in self.ddict:
                            self.ddict[sequence] = {}
                        if hmmprofile.startswith(self.leftprefix):
                            self._score(sequence, 'left', score, from_pos, to_pos, tlen)
                        elif hmmprofile.startswith(self.rightprefix):
                            self._score(sequence, 'right', score, from_pos, to_pos, tlen)
        except Exception as e:
            logging.error("Exception occurred when parsing HMMSearh results")
            raise e

    def __init__(self, domtable, region):
        self.domtable = domtable
        self.ddict = {}
        if region == "ITS2":
            self.leftprefix = '3_'
            self.rightprefix = '4_'
        elif region == "ITS1":
            self.leftprefix = '1_'
            self.rightprefix = '2_'
        elif region == "ALL":
            self.leftprefix = '1_'
            self.rightprefix = '4_'
        self.parse()


    def get_position(self, sequence):
        """ Returns the start and stop positions for a given sequence.
        Args:
            sequence (str): The name of the sequence.
        Returns:
            (tuple): (start position, end position) zero indexed
        Raises:
            KeyError: If input sequence is not present in dictionary (no ITS start or stop sites were found)
        """

        try:
            if "left" in self.ddict[sequence]:
                start = int(self.ddict[sequence]["left"]["to_pos"])
            else:
                start = None
            if "right" in self.ddict[sequence]:
                stop = int(self.ddict[sequence]["right"]["from_pos"]) - 1
            else:
                stop = None
            if "tlen" in self.ddict[sequence]:
                tlen = int(self.ddict[sequence]["tlen"])
            else:
                tlen = None
            return(start, stop, tlen)
        except KeyError:
            logging.debug("No ITS stop or start sites were identified for sequence {}, skipping.".format(sequence))
            raise KeyError
