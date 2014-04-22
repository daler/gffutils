from gffutils.helpers import asinterval

try:
    from pybedtools.contrib.plotting import Track
except ImportError:
    import warnings
    warnings.warn("Please install pybedtools for plotting.")


class Gene(object):
    def __init__(self, db, gene_id, transcripts=['mRNA'], utrs=["3'UTR",
        "5'UTR"], cds=['CDS'], ybase=0, **kwargs):
        """
        Represents a gene, `gene_id`, as a collection of
        pybedtools.contrib.plotting.Track objects.

        This class is flexible in how transcripts/CDSs/UTRs are defined;
        for example, you will usually have to adjust what featuretypes to
        consider as UTRs based on your annotation.

        For example, GTF files will typically define "3'UTR" but FlyBase GFF
        files use "three_prime_UTR".  In the former case you'd specify
        `utrs=["3'UTR", "5'UTR"]` while in the latter,
        `utrs=["three_prime_UTR", "five_prime_UTR"]`.

        The full length of the transcript is represented by a thin line; UTRs
        are slightly thicker, and CDSs are thickest.

        `kwargs` are passed to pybedtools.contrib.plotting.Track.

        `ybase` sets the bottom edge of the bottom isoform (sorted so that it's
        the shortest)

        Use the `max_y` attribute if you need to know the upper extent of the
        plotted transcripts.  This is useful if you are plotting multiple genes
        on the same axes and need to know what to use for the next one's
        `ybase`.

        You can adjust the `heights` attribute to set how wide transcripts,
        UTRs, CDSs are.  Padding is essentially "full" minus the largest height
        (CDS, 0.9, by default).
        """

        self.heights = {
                'transcript': 0.2,
                      'utrs': 0.5,
                       'cds': 0.9,
                      'full': 1.0}
        self.kwargs = kwargs
        self._transcripts = []
        for transcript in db.children(gene_id, level=1):
            if transcripts is None or transcript.featuretype in transcripts:
                d = {}
                d['transcript'] = [transcript]
                d['utrs'] = []
                d['cds'] = []
                for child in db.children(transcript, level=1):
                    _utrs = []
                    if child.featuretype in utrs:
                        d['utrs'].append(child)
                    if child.featuretype in cds:
                        d['cds'].append(child)
                self._transcripts.append(d)

        self.tracks = []

        self.ybase = ybase
        for d in sorted(self._transcripts,
                reverse=True, key=lambda x: len(x['transcript'])):
            self.tracks.append(self._make_track(d, 'transcript'))
            self.tracks.append(self._make_track(d, 'utrs'))
            self.tracks.append(self._make_track(d, 'cds'))
            self.ybase += self.heights['full']

        self.max_y = ybase + self.heights['full']

    def add_to_ax(self, ax):
        """
        Add all the transcripts for this gene to axes `ax`
        """
        for track in self.tracks:
            ax.add_collection(track)

    def _make_track(self, d, cls):
        yheight = self.heights[cls]
        ybase = self.ybase + (self.heights['full'] - yheight) * 0.5
        return Track(
                (asinterval(i) for i in d[cls]),
                ybase=ybase, yheight=yheight, **self.kwargs)
