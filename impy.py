import click
from path import Path
import subprocess
import os
import getpass
import json
import shlex
import urllib.request
import shutil
from copy import deepcopy
import tempfile
import sys
import re

CURRENT_DIRECTORY = Path(os.getcwd())
IMP_VERSION = '1.4'
IMP_DEFAULT_TAR_REPOSITORY = 'https://webdav-r3lab.uni.lu/public/R3lab/IMP/dist/imp-%s.tar.gz' % IMP_VERSION
IMP_IMAGE_NAME = 'docker-r3lab.uni.lu/imp/imp'
IMP_DEFAULT_DB_DIR = CURRENT_DIRECTORY / 'imp-db'
IMP_DEFAULT_OUTPUT_DIR = CURRENT_DIRECTORY / 'imp-output'
IMP_DEFAULT_CONFIG_FILE = CURRENT_DIRECTORY / 'userconfig.imp.json'

CONTAINER_OUTPUT_DIR = '/home/imp/output'
CONTAINER_DB_DIR = '/home/imp/databases'
CONTAINER_DATA_DIR = '/home/imp/data'
CONTAINER_CONF_DIR = '/home/imp/conf'
CONTAINER_CODE_DIR = '/home/imp/code'

"""
from functools import update_wrapper

def pass_obj(f):
    @click.pass_context
    def new_func(ctx, *args, **kwargs):
        return ctx.invoke(f, ctx.obj, *args, **kwargs)
    return update_wrapper(new_func, f)
"""

@click.group(context_settings={
    'help_option_names': ['-h', '--help']
})
@click.option('--image-name', default=IMP_IMAGE_NAME, help='IMP image name')
@click.option('--image-tag', default=IMP_VERSION, help='IMP image tag/version')
@click.option('--image-repo', default=IMP_DEFAULT_TAR_REPOSITORY, help='IMP repository. Must point to a gzipped image.')
@click.option('--enter', default=False, is_flag=True, help='Enter IMP docker image.')
@click.option('-d', '--database-path', help='Set different database path.', default=IMP_DEFAULT_DB_DIR)
@click.option('-c', '--config-file-path', help='Set different config file path.', default=IMP_DEFAULT_CONFIG_FILE)
@click.option('-s', '--source-code', help='Use IMP source code at the file path specified instead of the one shipped inside the image.')
@click.pass_context
def cli(ctx, image_name, image_tag, image_repo, database_path, config_file_path, source_code, enter):
    """Integrated Metaomic Pipeline"""
    if not ctx.obj:
        ctx.obj = {}
    # set shared context
    ctx.obj['image-name'] = image_name
    ctx.obj['image-tag'] = image_tag
    ctx.obj['image-repo'] = image_repo
    ctx.obj['database-path'] = Path(database_path).abspath()
    ctx.obj['config-file-path'] = Path(config_file_path).abspath()
    ctx.obj['source-code'] = Path(source_code).abspath()
    ctx.obj['enter'] = enter
    # validate
    if source_code is not None:
        source_code = Path(source_code)
        if not source_code.isdir():
            click.secho("`source code` must be a directory.", fg='red', bold=True)
            ctx.abort()
    database_path = Path(database_path)
    if not database_path.exists():
        database_path.makedirs()
    if not database_path.isdir():
        click.secho("`database path` must be a directory.", fg='red', bold=True)
        ctx.abort()

    # # TODO # add environment variables
    # envs = ['-e {}="{}"'.format(*e.split('=')) for e in args['-e']]

    # # TODO # split snakemake workflows into multiple

    # # TODO # workflow for single omics
@cli.command()
@click.option('--single-omics', is_flag=True, default=False, help='Activate single omics mode.')
@click.pass_context
def test(ctx, single_omics, *args):
    print(ctx.obj)
    print(single_omics)

def requirements():
    """
    Check if requirements are installed.
    """
    good = True
    # docker
    try:
        subprocess.check_output(['which', 'docker'])
    except subprocess.CalledProcessError:
        good = False
        click.secho("Docker must be installed. Please see https://docs.docker.com/installation.", fg='red', bold=True)
    try:
        subprocess.check_output(['docker', 'ps'])
    except subprocess.CalledProcessError:
        good = False
        click.secho("Docker must be used without sudo. Please see https://docs.docker.com/engine/installation/linux/ubuntulinux/#/create-a-docker-group.", fg='red', bold=True)
    # python3
    if sys.version_info < (3, 0, 0):
        good = False
        click.secho("Python 3 or later must be installed. Please see https://www.python.org/downloads.", fg='red', bold=True)
    return good

def call(cmd, container_name):
    click.secho("""Executing IMP command:
    %s
    """ % cmd, fg='green')
    try:
        p = subprocess.Popen(cmd, shell=True)
        # Poll process for new output until finished
        # while True:
        #     nextline = p.stdout.readline()
        #     if not nextline and p.poll() is not None:
        #         break
        #     click.secho(str(nextline, 'utf-8'))
            # sys.stdout.write(nextline)
            # sys.stdout.flush()
        p.wait()
        sys.exit(p.returncode)
    except KeyboardInterrupt:
        click.secho('Keyboard interruption. Killing container...', fg='green')
        k = subprocess.Popen("docker rm -f %s" % container_name, shell=True)
        click.secho('killed.', fg='green')
        return p.terminate()

def is_imp_container_installed(name, tag):
    """
    Check if IMP is installed
    """
    if not requirements():
        click.secho("Abort.", fg='red', bold=True)
        return
    try:
        imp_images = subprocess.Popen(['docker', 'images', name], stdout=subprocess.PIPE)
        is_installed = subprocess.check_output(['grep', " %s " % tag], stdin=imp_images.stdout)
    except subprocess.CalledProcessError:
        return False
    click.secho("[x] Found IMP {name} {tag}".format(name=name, tag=tag), fg='green')
    return True


@cli.command()
@click.pass_context
def install_imp_container(ctx):
    """
    Install IMP container.
    """
    if not requirements():
        click.secho("Abort.", fg='red', bold=True)
        return

    if is_imp_container_installed(ctx.obj['image-name'], ctx.obj['image-tag']):
        click.secho("Abort.", fg='red', bold=True)
        return
    fname = 'imp-tarball.tmp.tgz'
    # get image
    if  ctx.obj['image-repo'][:4].startswith('http'):  # download
        click.secho("[x] Downloading IMP TARBALL at '%s'" %  ctx.obj['image-repo'], fg='green')
        with urllib.request.urlopen(ctx.obj['image-repo']) as response, open(fname, 'wb') as out_file:
            shutil.copyfileobj(response, out_file)
    else: # copy
        click.secho("[x] Copying IMP TARBALL '%s'" %  ctx.obj['image-repo'], fg='green')
        with open(fname, 'wb') as out_file:
            shutil.copy(ctx.obj['image-repo'], out_file)
    try:
        # load image
        click.secho("[x] Loading IMP TARBALL into docker", fg='green')
        subprocess.check_output(['docker', 'load', '-i', fname])
    except:
        pass
    finally:
        # clean
        click.secho("[x] Removing IMP TARBALL.", fg='green')
        os.remove(fname)

def generate_container_name(directory):
    """
    Slugify container name base on the directory given
    """
    slug = re.sub('\W', '_', directory)
    while slug.startswith('_'):
        slug = slug[1:]
    return slug

def generate_docker_cmd(container_name, database_path, configuration_file_path,
                        image_name, image_tag, interactive,
                        data_directory=None,
                        command=None, source_code=None,
                        output_directory=None, environment=None):

    configuration_file_path = Path(configuration_file_path)
    configuration_file_name = str(Path(configuration_file_path).name)
    configuration_file_dir = Path(configuration_file_path).parent.abspath()

    # prepare general command
    cmd = "docker run --rm --name {container_name}".format(container_name=container_name)

    # add volumes
    volumes = " -v {database_path}:{container_database_dir} -v {configuration_file_dir}:{container_configuration_file_dir}".format(
        database_path=database_path,
        container_database_dir=CONTAINER_DB_DIR,
        container_configuration_file_dir=CONTAINER_CONF_DIR,
        configuration_file_dir=configuration_file_dir,
    )
    if source_code is not None:
        volumes += " -v {source_code}:{container_source_code_dir}".format(
            source_code=source_code,
            container_source_code_dir=CONTAINER_CODE_DIR
        )
    if output_directory is not None:
        volumes += " -v {output_directory}:{container_output_dir}".format(
            output_directory=output_directory,
            container_output_dir=CONTAINER_OUTPUT_DIR
        )
    if data_directory is not None:
        volumes += " -v {data_directory}:{container_data_dir}".format(
            data_directory=data_directory,
            container_data_dir=CONTAINER_DATA_DIR
        )
    cmd += volumes

    # add environment variables
    environments = """ -e "LOCAL_USER_ID=`id -u $USER`" -e "LOCAL_GROUP_ID=`id -g $USER`" -e CONFIGFILE={container_configuration_file_dir}/{configuration_file_name}""".format(
        container_database_dir=CONTAINER_DB_DIR,
        container_configuration_file_dir=CONTAINER_CONF_DIR,
        configuration_file_path=configuration_file_path,
        configuration_file_name=configuration_file_name
    )
    if environment is not None:
        for k, v in environment.items():
            environments += """ -e {key}="{value}" """.format(
                key=k,
                value=v
            )
    cmd += environments

    # add it flag if specified
    if interactive:
        cmd += " -it"
    # add container name:tag
    cmd += " {image_name}:{image_tag}".format(
        image_name=image_name,
        image_tag=image_tag
    )
    # add command
    if interactive:
        command = '/bin/bash'
    cmd += ' %s' % command
    return cmd


@cli.command()
@click.pass_context
def init(ctx):
    """
    Initialise databases
    """
    if not is_imp_container_installed(ctx.obj['image-name'], ctx.obj['image-tag']):
        click.secho('IMP image not installed. Please run `impy install_imp_container` first.', bold=True)

    container_name = generate_container_name(ctx.obj['database-path'])

    init_cmd = "snakemake -s {container_source_code_dir}/rules/ini/init".format(
        container_source_code_dir=CONTAINER_CODE_DIR
    )

    docker_cmd = generate_docker_cmd(
        container_name,
        ctx.obj['database-path'],
        ctx.obj['config-file-path'],
        image_name=ctx.obj['image-name'],
        image_tag=ctx.obj['image-tag'],
        command=init_cmd,
        interactive=ctx.obj['enter'],
        source_code=ctx.obj['source-code'],
        )

    call(docker_cmd, container_name)





@cli.command()
@click.option('-m', '--metagenomic', help="Path to the Metagenomic files.", multiple=True)
@click.option('-t', '--metranscriptomic', help="Path to the Metatranscriptomic files.", multiple=True)
@click.option('-o', '--output-directory', help="Output directory.", default=IMP_DEFAULT_OUTPUT_DIR)
@click.option('-a', '--assembler', help="Assembler to use.",  type=click.Choice(['idba', 'megahit']), default='megahit')
@click.option('--single-omics', is_flag=True, default=False, help='Activate single omics mode.')
@click.option('-x', '--execute',
              help="Command to execute.",
              default="snakemake -s {container_source_code_dir}/Snakefile".format(
              container_source_code_dir=CONTAINER_CODE_DIR))
@click.pass_context
def run(ctx, metagenomic, metranscriptomic,
        assembler, output_directory, single_omics,
        execute):
    """
    Run IMP workflow.

    Preprocessing --> Assembly --> Analysis --> Binning --> Report
    """

    # database path
    database_path = Path(ctx.obj['database-path']).abspath()

    # environment variable
    steps = ['preprocessing', 'assembly', 'analysis', 'binning', 'report']

    data_directory = None

    # find minimum common path between the data files
    # inorder to mount them in the container
    mg_data = [Path(p).abspath() for p in metagenomic]
    mt_data = [Path(p).abspath() for p in metranscriptomic]
    # check if paths exists
    for pth in mg_data + mt_data:
        if not pth.exists():
            click.secho('Path provided does not exists: `%s`.' % pth, fg='red', bold=True)
            ctx.abort()
    common_path = Path(os.path.commonprefix(mg_data + mt_data)).dirname()

    # update data paths to remove the 'common path' from it.
    mg_data = [p.partition(common_path)[-1][1:] for p in mg_data]
    mt_data = [p.partition(common_path)[-1][1:] for p in mt_data]
    # update data path to put the container path before
    mg_data = [CONTAINER_DATA_DIR + '/' + d for d in mg_data]
    mt_data = [CONTAINER_DATA_DIR + '/' + d for d in mt_data]

    # <-- preprocessing validation
    if 'preprocessing' in steps:
        # validate data input
        if single_omics:
            if mg_data and mt_data:
                click.secho('In `single omics` you should only provide `metagenomics` or `metatranscriptomics` data.', fg='red', bold=True)
                ctx.abort()
            if not mg_data and not mt_data:
                click.secho('In `single omics` you should provide `metagenomics` or `metatranscriptomics` data.', fg='red', bold=True)
                ctx.abort()
        else:
            if not mg_data or not mt_data:
                click.secho('You should provide `metagenomics` and `metatranscriptomics` data.', fg='red', bold=True)
                ctx.abort()
        if mg_data and len(mg_data) != 2:
            click.secho('Metagenomic data should be 2 paired files.', fg='red', bold=True)
            ctx.abort()
        if mt_data and len(mt_data) != 2:
            click.secho('Metatranscriptomic data should be 2 paired files.', fg='red', bold=True)
            ctx.abort()
        data_directory = common_path
    # <-- end preprocessing validation

    # <-- assembly without preprocessing validation
    if 'assembly' in steps and steps.index('assembly') == 0:
        if single_omics:
            if mg_data and mt_data:
                click.secho('In `single omics` you should only provide `metagenomics` or `metatranscriptomics` data.', fg='red', bold=True)
                ctx.abort()
            if not mg_data and not mt_data:
                click.secho('In `single omics` you should provide `metagenomics` or `metatranscriptomics` data.', fg='red', bold=True)
                ctx.abort()
        else:
            if not mg_data or not mt_data:
                click.secho('You should provide `metagenomics` and `metatranscriptomics` data.', fg='red', bold=True)
                ctx.abort()
        if mg_data and len(mg_data) != 3:
            click.secho('Metagenomic data should be 2 paired files and single end', fg='red', bold=True)
            ctx.abort()
        if mt_data and len(mt_data) != 3:
            click.secho('Metatranscriptomic data should be 2 paired files and single end', fg='red', bold=True)
            ctx.abort()
        data_directory = common_path
    # <-- end assembly validation

    ev = {
        'MG': ' '.join(mg_data),
        'MT': ' '.join(mt_data),
        'IMP_ASSEMBLER': assembler,
        'IMP_STEPS': ' '.join(steps)
    }

    # output directory
    output_directory = Path(output_directory).abspath()

    if not output_directory.exists():
        output_directory.makedirs()
    if not output_directory.isdir():
        click.secho("`output directory` must be a directory.", fg='red', bold=True)
        ctx.abort()

    container_name = generate_container_name(output_directory)

    run_cmd = "snakemake -s {container_source_code_dir}/Snakefile".format(
        container_source_code_dir=CONTAINER_CODE_DIR
    )
    if execute:
        run_cmd = execute
    # docker command
    docker_cmd = generate_docker_cmd(
        container_name,
        ctx.obj['database-path'],
        ctx.obj['config-file-path'],
        data_directory=data_directory,
        image_name=ctx.obj['image-name'],
        image_tag=ctx.obj['image-tag'],
        interactive=ctx.obj['enter'],
        source_code=ctx.obj['source-code'],
        command=run_cmd,
        output_directory=output_directory,
        environment=ev
        )

    # execute the command
    call(docker_cmd, container_name)


@cli.command()
@click.option('-m', '--metagenomic', help="Path to the Metagenomic files.", multiple=True)
@click.option('-t', '--metranscriptomic', help="Path to the Metatranscriptomic files.", multiple=True)
@click.option('-o', '--output-directory', help="Output directory.", default=IMP_DEFAULT_OUTPUT_DIR)
@click.option('-a', '--assembler', help="Assembler to use.",  type=click.Choice(['idba', 'megahit']), default='megahit')
@click.option('--single-omics', is_flag=True, default=False, help='Activate single omics mode.')
@click.option('-x', '--execute',
              help="Command to execute.",
              default="snakemake -s {container_source_code_dir}/Snakefile".format(
              container_source_code_dir=CONTAINER_CODE_DIR))
@click.option('--single-step', help="Only execute preprocessing.", is_flag=True)
@click.pass_context
def preprocessing(ctx, metagenomic, metranscriptomic,
        assembler, output_directory, single_omics,
        execute, single_step):
    """
    Run IMP workflow.

    Preprocessing --> Assembly --> Analysis --> Binning --> Report
    """

    # database path
    database_path = Path(ctx.obj['database-path']).abspath()

    # environment variable
    steps = ['preprocessing', 'assembly', 'analysis', 'binning', 'report']
    if single_step:
        steps = ['preprocessing']

    data_directory = None

    # find minimum common path between the data files
    # inorder to mount them in the container
    mg_data = [Path(p).abspath() for p in metagenomic]
    mt_data = [Path(p).abspath() for p in metranscriptomic]
    # check if paths exists
    for pth in mg_data + mt_data:
        if not pth.exists():
            click.secho('Path provided does not exists: `%s`.' % pth, fg='red', bold=True)
            ctx.abort()
    common_path = Path(os.path.commonprefix(mg_data + mt_data)).dirname()

    # update data paths to remove the 'common path' from it.
    mg_data = [p.partition(common_path)[-1][1:] for p in mg_data]
    mt_data = [p.partition(common_path)[-1][1:] for p in mt_data]
    # update data path to put the container path before
    mg_data = [CONTAINER_DATA_DIR + '/' + d for d in mg_data]
    mt_data = [CONTAINER_DATA_DIR + '/' + d for d in mt_data]

    # validate data input
    if single_omics:
        if mg_data and mt_data:
            click.secho('In `single omics` you should only provide `metagenomics` or `metatranscriptomics` data.', fg='red', bold=True)
            ctx.abort()
        if not mg_data and not mt_data:
            click.secho('In `single omics` you should provide `metagenomics` or `metatranscriptomics` data.', fg='red', bold=True)
            ctx.abort()
    else:
        if not mg_data or not mt_data:
            click.secho('You should provide `metagenomics` and `metatranscriptomics` data.', fg='red', bold=True)
            ctx.abort()
    if mg_data and len(mg_data) != 2:
        click.secho('Metagenomic data should be 2 paired files.', fg='red', bold=True)
        ctx.abort()
    if mt_data and len(mt_data) != 2:
        click.secho('Metatranscriptomic data should be 2 paired files.', fg='red', bold=True)
        ctx.abort()


    ev = {
        'MG': ' '.join(mg_data),
        'MT': ' '.join(mt_data),
        'IMP_ASSEMBLER': assembler,
        'IMP_STEPS': ' '.join(steps)
    }

    # output directory
    output_directory = Path(output_directory).abspath()

    if not output_directory.exists():
        output_directory.makedirs()
    if not output_directory.isdir():
        click.secho("`output directory` must be a directory.", fg='red', bold=True)
        ctx.abort()

    container_name = generate_container_name(output_directory)

    run_cmd = "snakemake -s {container_source_code_dir}/Snakefile".format(
        container_source_code_dir=CONTAINER_CODE_DIR
    )
    if execute:
        run_cmd = execute
    # docker command
    docker_cmd = generate_docker_cmd(
        container_name,
        ctx.obj['database-path'],
        ctx.obj['config-file-path'],
        data_directory=common_path,
        image_name=ctx.obj['image-name'],
        image_tag=ctx.obj['image-tag'],
        interactive=ctx.obj['enter'],
        source_code=ctx.obj['source-code'],
        command=run_cmd,
        output_directory=output_directory,
        environment=ev
        )

    # execute the command
    call(docker_cmd, container_name)




@cli.command()
@click.option('-m', '--metagenomic', help="Path to the Metagenomic files.", multiple=True)
@click.option('-t', '--metranscriptomic', help="Path to the Metatranscriptomic files.", multiple=True)
@click.option('-o', '--output-directory', help="Output directory.", default=IMP_DEFAULT_OUTPUT_DIR)
@click.option('-a', '--assembler', help="Assembler to use.",  type=click.Choice(['idba', 'megahit']), default='megahit')
@click.option('--single-omics', is_flag=True, default=False, help='Activate single omics mode.')
@click.option('-x', '--execute',
              help="Command to execute.",
              default="snakemake -s {container_source_code_dir}/Snakefile".format(
              container_source_code_dir=CONTAINER_CODE_DIR))
@click.option('--single-step', help="Only execute assembly step.", is_flag=True)
@click.pass_context
def assembly(ctx, metagenomic, metranscriptomic,
        assembler, output_directory, single_omics,
        execute, single_step):
    """
    Run IMP workflow.

    Assembly --> Analysis --> Binning --> Report
    """

    # database path
    database_path = Path(ctx.obj['database-path']).abspath()

    # environment variable
    steps = ['assembly', 'analysis', 'binning', 'report']
    if single_step:
        steps = ['assembly']

    data_directory = None

    # find minimum common path between the data files
    # inorder to mount them in the container
    mg_data = [Path(p).abspath() for p in metagenomic]
    mt_data = [Path(p).abspath() for p in metranscriptomic]
    # check if paths exists
    for pth in mg_data + mt_data:
        if not pth.exists():
            click.secho('Path provided does not exists: `%s`.' % pth, fg='red', bold=True)
            ctx.abort()
    common_path = Path(os.path.commonprefix(mg_data + mt_data)).dirname()

    # update data paths to remove the 'common path' from it.
    mg_data = [p.partition(common_path)[-1][1:] for p in mg_data]
    mt_data = [p.partition(common_path)[-1][1:] for p in mt_data]
    # update data path to put the container path before
    mg_data = [CONTAINER_DATA_DIR + '/' + d for d in mg_data]
    mt_data = [CONTAINER_DATA_DIR + '/' + d for d in mt_data]


    if single_omics:
        if mg_data and mt_data:
            click.secho('In `single omics` you should only provide `metagenomics` or `metatranscriptomics` data.', fg='red', bold=True)
            ctx.abort()
        if not mg_data and not mt_data:
            click.secho('In `single omics` you should provide `metagenomics` or `metatranscriptomics` data.', fg='red', bold=True)
            ctx.abort()
    else:
        if not mg_data or not mt_data:
            click.secho('You should provide `metagenomics` and `metatranscriptomics` data.', fg='red', bold=True)
            ctx.abort()
    if mg_data and len(mg_data) != 3:
        click.secho('Metagenomic data should be 2 paired files and single end', fg='red', bold=True)
        ctx.abort()
    if mt_data and len(mt_data) != 3:
        click.secho('Metatranscriptomic data should be 2 paired files and single end', fg='red', bold=True)
        ctx.abort()

    # <-- end assembly validation

    ev = {
        'MG': ' '.join(mg_data),
        'MT': ' '.join(mt_data),
        'IMP_ASSEMBLER': assembler,
        'IMP_STEPS': ' '.join(steps)
    }

    # output directory
    output_directory = Path(output_directory).abspath()

    if not output_directory.exists():
        output_directory.makedirs()
    if not output_directory.isdir():
        click.secho("`output directory` must be a directory.", fg='red', bold=True)
        ctx.abort()

    container_name = generate_container_name(output_directory)

    run_cmd = "snakemake -s {container_source_code_dir}/Snakefile".format(
        container_source_code_dir=CONTAINER_CODE_DIR
    )
    if execute:
        run_cmd = execute
    # docker command
    docker_cmd = generate_docker_cmd(
        container_name,
        ctx.obj['database-path'],
        ctx.obj['config-file-path'],
        data_directory=common_path,
        image_name=ctx.obj['image-name'],
        image_tag=ctx.obj['image-tag'],
        interactive=ctx.obj['enter'],
        source_code=ctx.obj['source-code'],
        command=run_cmd,
        output_directory=output_directory,
        environment=ev
        )

    # execute the command
    call(docker_cmd, container_name)


@cli.command()
@click.option('--data-dir', help="Path to the data directory containing output files from previous IMP step.")
@click.option('-o', '--output-directory', help="Output directory.", default=IMP_DEFAULT_OUTPUT_DIR)
@click.option('--single-omics', is_flag=True, default=False, help='Activate single omics mode.')
@click.option('-x', '--execute',
              help="Command to execute.",
              default="snakemake -s {container_source_code_dir}/Snakefile".format(
              container_source_code_dir=CONTAINER_CODE_DIR))
@click.option('--single-step', help="Only execute analysis step.", is_flag=True)
@click.pass_context
def analysis(ctx, data_dir, output_directory, single_omics,
             execute, single_step):
    """
    Run IMP workflow.

    Analysis --> Binning --> Report
    """

    # database path
    database_path = Path(ctx.obj['database-path']).abspath()

    # environment variable
    steps = ['analysis', 'binning', 'report']
    if single_step:
        steps = ['analysis']


    data_dir = Path(data_dir).abspath()
    if not data_dir.isdir():
        click.secho('--data-dir must be a directtory.', fg='red', bold=True)
        ctx.abort()

    # look for mg and mt data
    # PREPROCESSING
    preprocessing_dir = data_dir.dirs('Preprocessing')
    if preprocessing_dir:
        preprocessing_dir = preprocessing_dir[0]
    else:
        click.secho("`Preprocessing directory` not present.", fg='red', bold=True)
        ctx.abort()
    mg_preprocessing_files = ('mg.r1.fq', 'mg.r1.preprocessed.fq', 'mg.r2.fq', 'mg.r2.preprocessed.fq', 'mg.se.preprocessed.fq')
    mt_preprocessing_files = ('mt.r1.fq', 'mt.r1.preprocessed.fq','mt.r2.fq', 'mt.r2.preprocessed.fq', 'mt.se.preprocessed.fq')

    got_mg = True
    got_mt = True
    for f in mg_preprocessing_files:
        if not preprocessing_dir.files(f):
            got_mg = False
            break
            # click.secho("`Preprocessing directory` must contains '%s'." % f, fg='red', bold=True)
            # ctx.abort()
    for f in mt_preprocessing_files:
        if not preprocessing_dir.files(f):
            got_mt = False
            break
            # click.secho("`Preprocessing directory` must contains '%s'." % f, fg='red', bold=True)
            # ctx.abort()
    if single_omics:
        if not got_mg and not got_mt:
            click.secho("`Preprocessing directory` must contains mg or mt data. e.g: %s" % (', '.join(mg_preprocessing_files)), fg='red', bold=True)
            ctx.abort()
    else:
        if not got_mg or not got_mt:
            click.secho("`Preprocessing directory` must contains mg and mt data e.g: %s" % (', '.join(mg_preprocessing_files + mt_preprocessing_files)), fg='red', bold=True)
            ctx.abort()
    # update data path
    mg_data = [CONTAINER_DATA_DIR + '/' + d for d in mg_preprocessing_files]
    mt_data = [CONTAINER_DATA_DIR + '/' + d for d in mt_preprocessing_files]


    # ASSEMBLY
    assembly_dir = data_dir.dirs('Assembly')
    if assembly_dir:
        assembly_dir = assembly_dir[0]
    else:
        click.secho("`Assembly directory` not present.", fg='red', bold=True)
        ctx.abort()
    if single_omics:
        if not (assembly_dir.files('mg.assembly.merged.fa') and assembly_dir.files('mg.reads.sorted.bam')):
            if not (assembly_dir.files('mt.assembly.merged.fa') and assembly_dir.files('mt.reads.sorted.bam')):
                click.secho("`Assembly directory` must contains mg or mt data. e.g: mg.assembly.merged.fa, mg.reads.sorted.bam", fg='red', bold=True)
                ctx.abort()
    elif not (assembly_dir.files('mgmt.assembly.merged.fa')
              and assembly_dir.files('mt.assembly.merged.fa')
              and assembly_dir.files('mt.reads.sorted.bam')):
        click.secho("`Assembly directory` must contains mg and mt data. e.g: mgmt.assembly.merged.fa, mg.reads.sorted.bam, mg.reads.sorted.bam", fg='red', bold=True)
        ctx.abort()

    ev = {
        'MG': ' '.join(mg_data),
        'MT': ' '.join(mt_data),
        'IMP_STEPS': ' '.join(steps)
    }

    # output directory
    output_directory = Path(output_directory).abspath()

    if not output_directory.exists():
        output_directory.makedirs()
    if not output_directory.isdir():
        click.secho("`output directory` must be a directory.", fg='red', bold=True)
        ctx.abort()

    container_name = generate_container_name(output_directory)

    run_cmd = "snakemake -s {container_source_code_dir}/Snakefile".format(
        container_source_code_dir=CONTAINER_CODE_DIR
    )
    if execute:
        run_cmd = execute
    # docker command
    docker_cmd = generate_docker_cmd(
        container_name,
        ctx.obj['database-path'],
        ctx.obj['config-file-path'],
        data_directory=data_dir,
        image_name=ctx.obj['image-name'],
        image_tag=ctx.obj['image-tag'],
        interactive=ctx.obj['enter'],
        source_code=ctx.obj['source-code'],
        command=run_cmd,
        output_directory=output_directory,
        environment=ev
        )

    # execute the command
    call(docker_cmd, container_name)


if __name__ == '__main__':
    cli(obj={})
