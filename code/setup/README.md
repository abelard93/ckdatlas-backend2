## RUNNING setup.R

1. change wd_path to /your/path/to/git_repo/data/
2. change log_path to path in your directory where you want the log file to be written to 
3. make sure you have your set-url to git@gitlab.com:sysdiab/neo4j.git (not https)
    - this ensures that you are not prompted for your username/password 
4. make sure db_addr is correct
5. if the neo4j instance has been moved to a different server or you are running the script for the first time/from a different account:
    - update neo4j_server (ssh command to connect to server)
        - passwordless connection!
    - update wd_path_server /server/path/to/git_repo/


## LOG FILE

time stamp 	\t	node type		\t 	code 	\t 	free comment

code:

- FORMAT
- MAPPING
- VALUE
- INFO
