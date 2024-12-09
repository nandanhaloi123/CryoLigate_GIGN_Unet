import os
import sys
from datetime import datetime, timezone
from chimerax.core.commands import run


def log(
    message,
    status="ERROR",
    log_path=os.getcwd() + os.path.sep + "logs",
    log_filename="log.txt",
):
    """
    Creates logs.

    Args:
        message - message to put in the log file
        status - status of the message. Two options: INFO and ERROR
        log_path - path to the directory with log files
        log_filename - name of the log file
    """

    # full path to the log file (including its name)
    path_full = log_path + os.path.sep + log_filename

    # convert message to string in the case if exception was provided
    message = str(message)

    if not message.endswith("\n"):
        message += "\n"

    # log_message includes utc time of the message, status of the message
    # and the message itself
    log_message = (
        datetime.now(timezone.utc).strftime("%H:%M:%S.%f")
        + " "
        + status
        + ": "
        + message
    )

    with open(path_full, "a") as f:
        f.write(log_message)


def run_chimeraX_command_without_log(session, command):
    """
    Runs command inside ChimeraX without logging.

    Args:
        session - current ChimeraX session, available as a global variable when run
        python scripts inside ChimeraX
        command - string representing a command for ChimeraX
    """

    run(session, command)


def run_chimeraX_command_with_log(
    session,
    command,
    log_path=os.path.join(os.getcwd(), "chimeraX_logs"),
    log_filename="log.txt",
):
    """
    Runs command inside Chimera with logging.

    Args:
        session - current ChimeraX session, available as a global variable when run
        python scripts inside ChimeraX
        command - string representing a command for ChimeraX
        log_path - path to the directory with log files
        log_filename - name of the log file
    """

    try:
        run_chimeraX_command_without_log(session, command)
        log(
            f"Successfully executed command: {command}",
            status="INFO",
            log_path=log_path,
            log_filename=log_filename,
        )

    except Exception as e:
        log(
            f"Failed to execute the command: {command}. {e}",
            status="ERROR",
            log_path=log_path,
            log_filename=log_filename,
        )
        raise Exception(e)


def run_chimeraX_command(
    session,
    command,
    is_log=True,
    log_path=os.path.join(os.getcwd(), "chimeraX_logs"),
    log_filename="log.txt",
):
    """
    Runs command inside Chimera and logs results if needed.

    Args:
        session - current ChimeraX session, available as a global variable when run
        python scripts inside ChimeraX
        command - string representing a command for ChimeraX
        is_log - whether we write logs for commands executing
        log_path - path to the directory with log files
        log_filename - name of the log file
    """

    if is_log:
        run_chimeraX_command_with_log(
            session, command, log_path=log_path, log_filename=log_filename
        )

    else:
        run_chimeraX_command_without_log(session, command)


def delete_extension_from_filename(filename):
    """
    Deletes extension (symbols after .) from the given filename.

    Args:
        filename - name of the given file

    Returns:
        name of the file without extension
    """

    return ".".join(filename.split(".")[:-1])


def extract_filename_from_full_path(path_full):
    """
    Extracts filename from the full (including the name) path to file.

    Args:
        full_path - full path to the file

    Returns:
        filename extracted from the full path
    """

    return path_full.split(os.path.sep)[-1]
