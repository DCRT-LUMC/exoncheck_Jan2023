Flaskr
======

ExonCheck is built on top of the Flask `tutorial`_.

.. _tutorial: https://flask.palletsprojects.com/tutorial/

Important notes
-------
Domain information is not provided correctly yet, it shows general domain information and not the specific domains of the exon to be skipped
Currently working on showing the number of hits in the LOVD database

Install
-------
Install::

    # clone the repository
    $ git clone https://github.com/pallets/flask
    $ cd flask
    # checkout the correct version
    $ git tag  # shows the tagged versions
    $ git checkout latest-tag-found-above
    $ cd examples/tutorial

Create a virtualenv and activate it::

    $ python3 -m venv venv
    $ . venv/bin/activate

Or on Windows cmd
Download zip file, unzip file and go to directory in Windows Powershell. Subsequently do::

    > py -3 -m venv venv
    > venv\Scripts\activate.bat


Install Flaskr::

    > py -m pip install -e .

Or if you are using the main branch, install Flask from source before
installing Flaskr::

    > py -m pip install -e ../..
    > py -m pip install -e .

Install other required packages if not installed already::

    > py -m pip install -r requirements.txt


Run
---

::

Initialiase the database (only the first time or if you want a new database, **IT OVERWRITES THE EXISTING DATABASE**)::

    $ flask init-db

Run::

    $ export FLASK_APP=flaskr
    $ export FLASK_ENV=development
    $ flask run

Or on Windows cmd::

Initialiase the database (only the first time or if you want a new database, **IT OVERWRITES THE EXISTING DATABASE**)::

    > py -m flask --app flaskr init-db

Run::

    > py -m flask --app flaskr --debug run

Open http://127.0.0.1:5000 in a browser.


Test
----

::

    $ pip install '.[test]'
    $ pytest

Run with coverage report::

    $ coverage run -m pytest
    $ coverage report
    $ coverage html  # open htmlcov/index.html in a browser
