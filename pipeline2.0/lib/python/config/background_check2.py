#!/usr/bin/env python
import config_types

background = config_types.ConfigList('background')
background.add_config('screen_output', config_types.BoolConfig())
background.add_config('jobtracker_db', config_types.DatabaseConfig())
background.add_config('sleep', config_types.IntConfig())

if __name__=='__main__':
    import background2 as configs
    background2.populate_configs(configs.__dict__)
    background2.check_sanity()
