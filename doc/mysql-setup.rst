Local UCSC database
====================

If you wish to run a local version of the `UCSC Genome Database <https://genome.ucsc.edu/goldenPath/help/mysql.html>`_,
for example because queries to external databases are not allowed in your network, or because you want to use a
specific version of the database, you can do so by following this guide.

First of all, to be able to run a local version of this database, MySQL is required on your machine. Install it by
following the steps in `https://dev.mysql.com/downloads/installer/ <https://dev.mysql.com/downloads/installer/>`_.

Download database tables
-------------------------

Once MySQL is installed, you need to download the UCSC database tables :code:`refGene` and :code:`snp150Common` to your
local machine. You can do so by running the following commands:

.. code-block:: bash

    mysqldump -v --single-transaction -u genome -h genome-euro-mysql.soe.ucsc.edu hg38 snp150Common > snp150Common.sql

.. code-block:: bash

    mysqldump -v --single-transaction -u genome -h genome-euro-mysql.soe.ucsc.edu hg38 refGene > refGene.sql

* :code:`-v` is to activate verbose mode (not strictly necessary)
* :code:`--single-transaction` is to avoid locking the tables
* :code:`-u genome` is to specify the user (genome)
* :code:`-h genome-euro-mysql.soe.ucsc.edu` is to specify the host
* :code:`hg38` is the genome version
* :code:`snp150Common` and :code:`refGene` are the tables to download
* :code:`snp150Common.sql` and :code:`refGene.sql` are the output files

Create and populate local database
----------------------------------

Access MySQL CLI

.. code-block:: bash

    sudo mysql

Create a new database :code:`hg38`

.. code-block:: sql

    CREATE DATABASE hg38;

Review database

.. code-block:: sql

    SHOW CREATE DATABASE hg38;

Select the database

.. code-block:: sql

    USE hg38;

Copy the content of the downloaded files to the database

.. code-block:: bash

    SOURCE refGene.sql

.. code-block:: bash

    SOURCE snp150Common.sql


Create new user to access the database
--------------------------------------

Access MySQL

.. code-block:: bash

    sudo mysql

Create a new user :code:`genome` with password :code:`password`

.. code-block:: sql

    CREATE USER 'genome'@'localhost' IDENTIFIED BY 'password'

Grant privileges to the user :code:`genome`

.. code-block:: sql

    GRANT ALL PRIVILEGES ON hg38.* TO 'genome'@'localhost' WITH GRANT OPTION;

