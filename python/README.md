## Setup

To install all the requirements run:

```
pip3 install -r requirements.txt
make all
```

For macOS: If you are getting ‘cmath not found’ error from compiling, then you should consider downgrading xcode. There was a case with macos Sequoia 15.1.1 where xcode=16.1 led to ‘cmath not found’ and downgrading it to 16.0 fixed it.

## Acknowledgments
We thank Ilja Kuzborskij for helping porting the Matlab code to Python and Nicolò Campolongo for suggesting improvements.
