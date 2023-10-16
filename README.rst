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
Mac / Linux:
~~~~~~~~~~~~

.. code-block:: bash

  # clone the ExonCheck repository
  git clone https://github.com/DCRT-LUMC/exoncheck_Jan2023.git

  # Enter the folder
  cd exoncheck_Jan2023

  # Create a vritual environment
  python3 -m venv venv

  # Activate the virtual environment
  source venv/bin/activate

  # Install the requirements
  pip3 install -r requirements.txt 

  # Install ExonCheck
  pip3 install -e .

  # Set the correct environment variables
  export FLASK_APP=flaskr
  export FLASK_ENV=development

  # Create the database (only do this the first time, or create a backup first)
  cp instance/flaskr.sqlite{,.backup}
  flask init-db

  # Start the tool
  flask run


Open http://127.0.0.1:5000 in a browser.

Windows:
~~~~~~~~

Download zip file, unzip file and go to directory in Windows Powershell. Subsequently do:

    > py -3 -m venv venv
    > venv\Scripts\activate.bat

Install Flaskr::

    > py -m pip install -e .

Install other required packages if not installed already::

    > py -m pip install -r requirements.txt

Initialiase the database (only the first time or if you want a new database, **IT OVERWRITES THE EXISTING DATABASE**)::
> cp instance\flaskr.sqlite instance\flaskr.sqlite.backup
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
