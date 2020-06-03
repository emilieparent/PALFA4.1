import time_vs_phase
import freq_vs_phase

class ProfileClass(time_vs_phase.TimeVsPhaseClass, freq_vs_phase.FreqVsPhaseClass):
    data_key = "profile"

    def _compute_data(self, cand):
        """Create the dedispersed optimal profile for the candidate.

            Input:
                cand: A ratings2.0 Candidate object.

            Output:
                prof: The corresponding Profile object.
        """
        pfd = cand.get_from_cache('pfd')
        tvph = cand.get_from_cache('time_vs_phase')
        if pfd.fold_pow == 1.0:
            bestp = pfd.bary_p1
            bestpd = pfd.bary_p2
            bestpdd = pfd.bary_p3
        else:
            bestp = pfd.topo_p1
            bestpd = pfd.topo_p2
            bestpdd = pfd.topo_p3
        tvph.adjust_period(bestp, bestpd, bestpdd)
        return tvph.get_profile()
