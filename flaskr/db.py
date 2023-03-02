# sqlite3 is a built-in module in Python to interact with the SQLite database
import sqlite3
import click

from flask import current_app
from flask import g
from flask.cli import with_appcontext


def get_db():
    """Connect to the application's configured database. The connection
    is unique for each request and will be reused if this is called
    again in the same session.

    The function is called to create a connection to the database.
    """
    # g is a special object that stores the ongoing connection
    if "db" not in g:
        # A connection is established to the sqlite database.
        g.db = sqlite3.connect(
            current_app.config["DATABASE"], detect_types=sqlite3.PARSE_DECLTYPES
        )
        # Tells the connection to return rows that behave like dicts
        g.db.row_factory = sqlite3.Row

    return g.db


def close_db(e=None):
    """If this request connected to the database (i.e. g.db is set), close the
    connection.

    It is called after every request to close down the connection automatically.
    """
    db = g.pop("db", None)

    if db is not None:
        db.close()


def init_db():
    """Clear existing data and create new tables."""

    # Returns a database connection
    db = get_db()

    with current_app.open_resource("schema.sql") as f:
        db.executescript(f.read().decode("utf8"))

# Defines a command line command to initiate tables
@click.command("init-db")
@with_appcontext
def init_db_command():
    """Clear existing data and create new tables."""
    init_db()
    click.echo("Initialized the database.")


def init_app(app):
    """Register database functions with the Flask app. This is called by
    the application factory.
    """
    # Tells Flask to call close_db when cleaning up after returning the response
    app.teardown_appcontext(close_db)
    # Adds a new command that can be called with the Flask command
    app.cli.add_command(init_db_command)
