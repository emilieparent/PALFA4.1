#!/usr/bin/env python
import config_types

download = config_types.ConfigList('download')
download.add_config('api_service_url', config_types.StrConfig())
download.add_config('api_username', config_types.StrConfig())
download.add_config('api_password', config_types.StrConfig())
download.add_config('datadir', config_types.DirConfig())
download.add_config('space_to_use', config_types.IntOrLongConfig())
download.add_config('min_free_space', config_types.IntOrLongConfig())
download.add_config('numdownloads', config_types.IntConfig())
download.add_config('numrestored', config_types.IntConfig())
download.add_config('numretries', config_types.IntConfig())
download.add_config('ftp_host', config_types.StrConfig())
download.add_config('ftp_port', config_types.IntConfig())
download.add_config('ftp_username', config_types.StrConfig())
download.add_config('ftp_password', config_types.StrConfig())
download.add_config('request_timeout', config_types.IntConfig())
download.add_config('request_numbits', config_types.IntConfig())
download.add_config('request_datatype', config_types.StrConfig())
download.add_config('use_lftp', config_types.BoolConfig())

if __name__=='__main__':
    import download as configs
    download.populate_configs(configs.__dict__)
    download.check_sanity()