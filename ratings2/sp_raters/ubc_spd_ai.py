import numpy as np
import scipy as sp

import presto

import base
from sp_rating_classes import profile
from ratings2 import config_rater 

#### setup UBC AI
from ubc_AI.data import pfdreader
import cPickle
clfer = config_rater.spd_classifier
classifier = cPickle.load(open(clfer, 'rb'))
####

class ubc_spd_ai(base.BaseRater):
    short_name = "spd_AI"
    long_name = "UBC spd AI"
    description = "compute the prediction from the classifier " \
                  "based on spd files."
    version = 4
    
    rat_cls = profile.ProfileClass()

    def _compute_rating(self, cand):
        """Return a rating for the candidate. The rating value is 
            the prepfold sigma value.

            Input:
                cand: A Candidate object to rate.

            Output:
                value: The rating value.
        """
        pred = classifier.report_score([pfdreader(cand.spdfn)])
        return pred


Rater = ubc_spd_ai 
