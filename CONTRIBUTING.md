# Welcome!
Thank you for contributing to CDC's Open Source projects! If you have any
questions or doubts, don't be afraid to send them our way. We appreciate all
contributions, and we are looking forward to fostering an open, transparent, and
collaborative environment.

Before contributing, we encourage you to also read or [LICENSE](LICENSE),
[README](README.md), and
[code-of-conduct](CODE_OF_CONDUCT.md)
files. If you have any inquiries or questions not
answered by the content of this repository, feel free to [contact us](mailto:xlr@cdc.gov).

## Public Domain
This project is in the public domain within the United States, and copyright and
related rights in the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this project will be released under the CC0 dedication. By
submitting a pull request or pushing code directly to a branch, you are agreeing to comply with this waiver of
copyright interest.

## Privacy
All comments, messages, pull requests, and other submissions received through
CDC including this GitHub repository are subject to the [Presidential Records Act](https://www.archives.gov/about/laws/presidential-records.html)
and may be archived. Learn more at [https://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Using the issue tracker

The issue tracker is our preferred channel for [bug reports](#bug-reports), [features requests](#feature-requests)
and [submitting pull requests](#pull-requests).

> Security vulnerabilities should be emailed to xlr@cdc.gov and should not be reported on the issue tracker. We'll get back to you within five business days. Please be as detailed and specific as possible when emailing such reports.

## Bug reports

Guidelines for bug reports:

1. **Search prior issues** to determine if the problem has already been reported.

2. **Check if the issue has been fixed** by trying to reproduce it on the master or develop branches.

3. **Include the steps needed to reproduce the issue**

A good bug report is as detailed as possible. It should include information about your environment such as your browser version, runtime, and operating system. It should include setps needed to reproduce the problem, what you expected the outcome of those steps to be, and if applicable, screenshots. These details will help contributors in fixing potential bugs.

## Feature requests

We welcome feature requests. Before submitting a feature request, please consider whether the idea fits within the scope and intent of the project. The onus is on the requestor to make a compelling case on the merits of the feature.

## Pull requests
Pull requests (PRs) for fixes, improvements, and features are great and we strongly encourage them. Please ensure PRs are scoped to the project and avoid unrelated commits.

**Please ask first** before embarking on a significant pull request. You may spend a great deal of time on changes that the team doesn't want to include or that someone else was already working on.

PRs that fail to adhere to the coding conventions, styles, and other requirements (e.g. test coverage) may be asked to address those problems prior to a merge or rejected outright.

Adhering to the following process is the best way to get your work
included in the project:

1. [Fork](https://help.github.com/articles/fork-a-repo/) the project, clone your
   fork, and configure the remotes:

   ```bash
   # Clone your fork of the repo into the current directory
   git clone https://github.com/<your-username>/nndss-synthetic-data-generator.git
   # Navigate to the newly cloned directory
   cd nndss-synthetic-data-generator
   # Assign the original repo to a remote called "upstream"
   git remote add upstream https://github.com/cdcent/nndss-synthetic-data-generator.git
   ```

2. If you cloned a while ago, get the latest changes from upstream:

   ```bash
   git checkout master
   git pull upstream master
   ```

3. Create a new topic branch (off the main project development branch) to
   contain your feature, change, or fix:

   ```bash
   git checkout -b <topic-branch-name>
   ```

4. Commit your changes in logical chunks. Please adhere to these [git commit
   message guidelines](https://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html)
   or your code is unlikely be merged into the main project. Use Git's
   [interactive rebase](https://help.github.com/articles/about-git-rebase/)
   feature to tidy up your commits before making them public.

5. Locally merge (or rebase) the upstream development branch into your topic branch:

   ```bash
   git pull [--rebase] upstream master
   ```

6. Push your topic branch up to your fork:

   ```bash
   git push origin <topic-branch-name>
   ```

7. [Open a Pull Request](https://help.github.com/articles/using-pull-requests/)
    with a clear title and description.