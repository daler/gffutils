import inspect
import psutil


def tracer(text="", lookfor=".gtf"):
    offset = 0
    frames = inspect.stack()[::-1]
    loc = []
    pad = 0
    for frame in frames[1:-1]:
        frame = frame[0]
        try:
            name = frame.f_code.co_name
            line = frame.f_lineno
            if 'self' in frame.f_locals:
                cls = frame.f_locals['self'].__class__.__name__
                name = '.'.join([cls, name])
            loc.append('{pad}{name} <{line}>'.format(
                pad=' ' * pad, name=name, line=line))
        except AttributeError:
            continue
        pad += 1
    loc = '\n-> '.join(loc)
    proc = psutil.Process()
    opened = sum(1 for i in proc.open_files() if lookfor in i.path)
    print('[{opened}] {text}\n{loc}\n'.format(**locals()))
