import scipy as sp
import numpy as np
import presto
import sys

import base
from rating_classes import profile
import config_rater

#### setup UBC AI
from ubc_AI.data import pfdreader
import numpy as np
import cPickle
clfer = config_rater.pfd_classifier
classifier = cPickle.load(open(clfer, 'rb'))
####

class ubc_pfd_ai(base.BaseRater):
    short_name = "pfd_AI"
    long_name = "UBC pfd AI"
    description = "compute the prediction from the pulsar classifier " \
                  "based on pfd files."
    version = 5
    
    rat_cls = profile.ProfileClass()

    def _compute_rating(self, cand):
        """Return a rating for the candidate. The rating value is 
            the prepfold sigma value.

            Input:
                cand: A Candidate object to rate.

            Output:
                value: The rating value.
        """
        pfd_fn = cand.get_from_cache('pfd').pfd_filename

        pred = classifier.report_score([pfdreader(pfd_fn)])
        return pred


Rater = ubc_pfd_ai 
