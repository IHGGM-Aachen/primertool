Setup
==========

There are two ways to run the streamlit interface to this project. The first is to run the streamlit interface locally,
which is the easiest way to get started and recommended for development, because changes on the code will be effective
immediately (and do not require rebuilding an image).

The second option is to run the streamlit interface in a docker container. This is recommended for production, because
it is easier to deploy and scale.

If you cannot access the UCSC MySQL database from your local machine, follow the instructions for creating a
:doc:`mysql-setup`.

Streamlit
---------

To run the streamlit interface locally, you need to have an environment with the required dependencies installed. It is
recommended to do so in a dedicated environment (e.g. using conda or virtualenv). To create a virtualenv, you can run

.. code-block:: bash

    python3 -m venv primertool_env
    source primertool_env/bin/activate
    pip install -r requirements.txt

or create a conda environment by running

.. code-block:: bash

    conda env create -f environment.yml

Then, you can simply start the streamlit interface by running

.. code-block:: bash

    streamlit run streamlit_main.py --server.port=8501 --server.address=0.0.0.0

in your environment. This will start a local server at :code:`http://localhost:8501`/:code:`http://0.0.0.0:8501` and
open a browser window with the interface.


Docker
------

To run the streamlit interface in a docker container, you need to have docker installed on your machine
(`install docker <https://docs.docker.com/engine/install/>`_). You can then build the docker image by running

.. code-block:: bash

    sudo docker build -t primertool .

in the root directory of this project. This will create a docker image called :code:`primertool`.

.. note::
    If you are behind a firewall and can't build the image because of a proxy error, you can set the proxy in the
    :code:`Dockerfile`. Simply uncomment the lines in the Dockerfile that set the proxy and set the proxy to your own proxy.

After briefly checking if the image was build with :code:`sudo docker images`, you can then run the image in a container
by running

.. code-block:: bash

    sudo docker run -p 8501:8501 primertool

or

.. code-block:: bash

    sudo docker run -it --net=host -p 8501:8501 primertool

if you are using a local version of UCSC to enable access to the local MySQL database.

To stop the container, you can run

.. code-block:: bash

    sudo docker stop <container_id>

where :code:`<container_id>` is the id of the container that you can find by running :code:`sudo docker ps`.

To start the container again, you can run

.. code-block:: bash

    sudo docker start <container_id>

Log access
----------

To gain more information on how the application is used and to monitor the application, activities as well as errors
are logged. If you run the application without docker, the logs are simply stored in the file
:code:`logs/primertool.log`. If you run the application in a docker container, the logs are stored in the container
and need to be copied onto your filesystem first. You can do so by running

.. code-block:: bash

    sudo docker cp <container_id>:/primertool.log /path/to/destination

where :code:`<container_id>` is the id of the container that you can find by running :code:`sudo docker ps` and
:code:`/path/to/destination` is the path on your filesystem where you want to store the log file.
