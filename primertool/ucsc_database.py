"""
This module handles the connection and queries to the UCSC SQL database.
"""
import mysql.connector
from mysql.connector import errorcode
import logging
from . import exceptions

logger = logging.getLogger(__package__)


def query(genome_assembly: str, query: str, local: bool = False, password_local: str = 'password') -> list or None:
    """Query the UCSC SQL database or a local copy of the database.

    Args:
        genome_assembly (str): The genome assembly to query, e.g. 'hg38'.
        query (str): The SQL query to execute.
        local (bool): If True, a local copy of the UCSC database is used.
        password_local (str): The password for the local database.
    Returns:
        list or None: The result of the query as a list of tuples or None if the query did not return any results.
    """
    # UCSC SQL database config
    ucsc_config = dict(user='genome',
                       password='',
                       host='genome-euro-mysql.soe.ucsc.edu',
                       database=genome_assembly,
                       raise_on_warnings=True,
                       )
    if local:
        ucsc_config['password'] = password_local
        ucsc_config['host'] = 'localhost'

    query_result = None
    try:
        with mysql.connector.connect(**ucsc_config) as connection, connection.cursor() as cursor:
            cursor.execute(query)
            query_result = cursor.fetchall()
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            logger.error("Database access denied. Check your username and password.")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            logger.error(f"Database '{genome_assembly}' does not exist.")
        else:
            logger.error(err)

    if query_result:
        return query_result
    else:
        logger.debug('This database query did not return any results. Please check your input.')
        return None
