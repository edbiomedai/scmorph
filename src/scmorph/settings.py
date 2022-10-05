settings = {}  # type: ignore
settings["cachedir"] = None


def get_cachedir() -> str:
    if not settings["cachedir"]:
        from tempfile import TemporaryDirectory

        settings["cachedir"] = TemporaryDirectory().name
    return settings["cachedir"]
