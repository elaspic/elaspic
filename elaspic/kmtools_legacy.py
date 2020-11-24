"""
Legacy stuff extracted from older versions of kmtools so that
we don't have to include kmtools as a dependency.
"""
import logging
import os
import os.path as op
import re
import shlex
import subprocess
from contextlib import contextmanager

logger = logging.getLogger(__name__)

# === System tools ===


@contextmanager
def decompress(filename):
    """Temporarly decompress a file."""
    try:
        logger.info("Gunzipping file '{}'...".format(filename))
        subprocess.check_call(shlex.split("gunzip '{}'".format(filename)))
    except Exception as e:
        logger.error("Unzipping the file failed with an error: {}".format(e))
        raise e
    else:
        yield op.splitext(filename)[0]
    finally:
        logger.info("Gzipping the file back again...")
        subprocess.check_call(shlex.split("gzip '{}'".format(filename.rstrip(".gz"))))


@contextmanager
def switch_paths(working_path):
    current_path = os.getcwd()
    try:
        os.chdir(working_path)
        yield
    except:
        raise
    finally:
        os.chdir(current_path)


# === DB tools ===


def parse_connection_string(connection_string):
    """Split `connection_string` into a dictionary of connection properties.

    .. note::

        The returned dictionary maps everything to strings (including `db_port`)
        to make it compatible with :class:`configparser.configParser`.

    Parameters
    ----------
    connection_string : str
        String describing database connection in SQLAlchemy-compatible format.

    Returns
    -------
    Mapping[str, str]

    Examples
    --------
    >>> from pprint import pprint
    >>> pprint(parse_connection_string('mysql://user:@localhost'))
    {'db_password': '',
     'db_port': '',
     'db_schema': '',
     'db_socket': '',
     'db_type': 'mysql',
     'db_url': 'localhost',
     'db_username': 'user'}
    >>> pprint(parse_connection_string('mysql://user:pass@192.168.0.1:3306/test'))
    {'db_password': 'pass',
     'db_port': '3306',
     'db_schema': 'test',
     'db_socket': '',
     'db_type': 'mysql',
     'db_url': '192.168.0.1',
     'db_username': 'user'}
    >>> pprint(parse_connection_string('sqlite:////absolute/path/to/foo.db'))
    {'db_password': None,
     'db_port': '',
     'db_schema': '/absolute/path/to/foo.db',
     'db_socket': '',
     'db_type': 'sqlite',
     'db_url': '',
     'db_username': ''}
    >>> connection_string = 'mysql://user@192.168.0.1:3306/test?unix_socket=/tmp/mysql.sock'
    >>> pprint(parse_connection_string(connection_string))
    {'db_password': None,
     'db_port': '3306',
     'db_schema': 'test',
     'db_socket': '/tmp/mysql.sock',
     'db_type': 'mysql',
     'db_url': '192.168.0.1',
     'db_username': 'user'}
    """
    db_params = {}
    (
        db_params["db_type"],
        db_params["db_username"],
        db_params["db_password"],
        db_params["db_url"],
        db_params["db_port"],
        db_params["db_schema"],
        db_params["db_socket"],
    ) = re.match(
        "^(\w*)"  # db_type
        "://"
        "(|\w*)"  # db_username
        "(|:\w*)"  # db_password
        "(|@localhost|@[a-zA-Z0-9\.-]*|@[0-9\.]*)"  # db_url
        "(|:[0-9]*)"  # db_port
        "(|\/[^?]*)"  # db_schema
        "(|\?unix_socket=.*)$",  # db_socket
        connection_string,
    ).groups()
    if db_params["db_password"].startswith(":"):
        if db_params["db_password"] == ":":
            db_params["db_password"] = ""
        else:
            db_params["db_password"] = db_params["db_password"][1:]
    else:
        db_params["db_password"] = None
    db_params["db_url"] = db_params["db_url"].lstrip("@")
    db_params["db_port"] = db_params["db_port"].lstrip(":")
    db_params["db_schema"] = (
        db_params["db_schema"][1:]
        if db_params["db_schema"].startswith("/")
        else db_params["db_schema"]
    )
    db_params["db_socket"] = db_params["db_socket"].partition("?unix_socket=")[-1]
    return db_params


def make_connection_string(**vargs):
    """Join a dictionary of connection properties (`vargs`) into a connection string.

    Examples
    --------
    >>> make_connection_string(**{ \
        'db_password': None, \
        'db_port': '', \
        'db_schema': '', \
        'db_socket': '', \
        'db_type': 'mysql', \
        'db_url': 'localhost', \
        'db_username': 'user'})
    'mysql://user@localhost/'
    >>> make_connection_string(**{ \
        'db_password': 'pass', \
        'db_port': '3306', \
        'db_schema': 'test', \
        'db_socket': '', \
        'db_type': 'mysql', \
        'db_url': '192.168.0.1', \
        'db_username': 'user'})
    'mysql://user:pass@192.168.0.1:3306/test'
    >>> make_connection_string(**{ \
        'db_password': '', \
        'db_port': '', \
        'db_schema': '/absolute/path/to/foo.db', \
        'db_socket': '', \
        'db_type': 'sqlite', \
        'db_url': '', \
        'db_username': ''})
    'sqlite:////absolute/path/to/foo.db'
    """
    vargs["db_password"] = (
        ":{}".format(vargs["db_password"])
        if vargs.get("db_password") is not None
        and not vargs.get("db_schema", "").startswith("/")
        else ""
    )
    vargs["db_url"] = "@{}".format(vargs["db_url"]) if vargs.get("db_url") else ""
    vargs["db_port"] = ":{}".format(vargs["db_port"]) if vargs.get("db_port") else ""
    vargs["db_schema"] = (
        "/{}".format(vargs["db_schema"]) if vargs.get("db_schema") else "/"
    )
    vargs["db_socket"] = (
        "?unix_socket={}".format(vargs["db_socket"]) if vargs.get("db_socket") else ""
    )
    connection_string = "{db_type}://{db_username}{db_password}{db_url}{db_port}{db_schema}{db_socket}".format(
        **vargs
    )
    return connection_string


@contextmanager
def lock_tables(tablenames, engine):
    """Lock a list of tables so that other processes can't access them while you do your thing."""
    if type(tablenames) not in {list, tuple}:
        tablenames = [tablenames]
    try:
        engine.execute("set innodb_lock_wait_timeout=14400")
        engine.execute("lock tables " + " ".join([t + " write" for t in tablenames]))
        yield
    except:
        raise
    finally:
        engine.execute("unlock tables")
