### Install requirements
`pip install -r requirements.txt`

### Running tests

Unit tests attempt to check that everything is working properly.
It is a good idea to run unit tests frequently, especially before making
changes and after making changes but before committing them.

run 
[nosetests](http://pythontesting.net/framework/nose/nose-introduction/) 
(this runs the unit tests):

```
nosetests
```

Output should look like this:

```
.......S...............
----------------------------------------------------------------------
Ran 23 tests in 2.179s

OK (SKIP=1)
```

Each dot stands for a unit test that ran, "S" stands for "Skipped".  If there are
failures the output will be more extensive, describing which tests failed and how.

For debugging purposes it is sometimes useful to use `print` statements and invoke
nosetests with the `--nocapture` option in order to see the output.
```
nosetests --nocapture
