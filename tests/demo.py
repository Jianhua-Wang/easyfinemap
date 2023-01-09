from rich.progress import Progress, TextColumn, BarColumn, MofNCompleteColumn, TimeElapsedColumn
from concurrent.futures import ProcessPoolExecutor
import time
from pathos.pools import _ProcessPool as Pool

def func(x, y, **kwargs):
    time.sleep(2)
    print(kwargs.get('a', 1))
    return x + y


args = [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6]]
args = [{'x': i, 'y': i+1, 'a': 2} for i in range(5)]

out = []
with Progress(
    TextColumn("{task.description}"),
    BarColumn(),
    MofNCompleteColumn(),
    TimeElapsedColumn(),
    auto_refresh=False,
) as progress:
    task = progress.add_task("Perform Fine-mapping...", total=len(args))
    with Pool(2) as p:
        results = [p.apply_async(func, kwds=kwargs) for kwargs in args]
        for result in results:
            progress.advance(task)
            out.append(result.get())
            progress.refresh()
print(out)
