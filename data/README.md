# Course Data Directory Instructions for Training Team
## Purpose
This directory contains all datasets and supplementary materials which are essential for the course. The files here support the learning objectives and provide the necessary resources for attendees.

## Suggested File Types
 - Datasets (CSV, BED+BIM+FAM, SAM, FASTA etc.)
 - Reference/supplementary data
 - Practical exercise data

## Best Practices for Naming Conventions
To ensure clarity and consistency, please follow these naming conventions for your files:
 - Use descriptive names that reflect the content (e.g., `dataset_experiment_1.csv`, `case-control.{bed,bim,fam}`).
 - Use lowercase letters and separate words with underscores (\_) for better readability.
 - Include version numbers if necessary (e.g., `data_v1.0.csv`).

## File Upload Guidelines
 - Files under 25MB can be uploaded via the web interface.
 - Files up to 100MB can be uploaded using the Git command line.
 - For files larger than 100MB, use Git Large File Storage (LFS) to manage your uploads.

## Steps to Add Data via GitHub Web Interface
1. Navigate to the Repository: Go to the GitHub repository for the course.
2. Select the **course_data** Directory: Click on the **course_data** folder to open it.
3. Click on "Add file": On the right-hand side, you'll see an option labelled Add file.
4. Choose "Upload files": From the dropdown menu, select Upload files.
5. Drag and Drop or Select Files: You can either drag and drop your files into the space provided or click choose your files to upload them from your computer.
6. Commit Changes: Once your files are uploaded, scroll down to the Commit changes section.
7. Enter a brief description of the changes in the commit message field (e.g., "Adding new dataset files").
8. Click Commit changes to finalize your upload.

## Git LFS Instructions
1. Install Git LFS.
2. Clone the repository.
3. Copy your large files to the respective **course_data** directory.
4. Track the suffix of large files using the following command:
```
git lfs track "*.sam"
```

5. Then add the .gitattributes, commit, and push:
```
git add .gitattributes
git commit -m "Now tracking sam files"
git push
```

6. Now add your files, then commit and push:
```
git add course_data/dir/file.sam
 git commit -m "Adding big sam file"
 git push
```

NB: Git LFS does NOT track changes in binary files--use different names for modified files.

## Informatics Software Guide
The Software Informatics Guide provides comprehensive information about the software utilized by participants throughout the course, including version numbers, links to software repositories or websites, and other relevant details. The training team is encouraged to review this guide and add their software, as this will be beneficial for course attendees in the future. The link to this guide can be found on the course homepage.

# Helpful Links
[Using Markdown on GitHub](https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github)
[Formatting markdown on GitHub](https://docs.github.com/en/github/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax)
[GitHub Markdown Cheatsheet](https://github.github.io/gfm/)
[Command Line GitHub Manual](https://cli.github.com/manual/)
[Git Large File Storage (LFS)](https://git-lfs.github.com/)

---

[African Genomics Short Course GitHub Home Page](https://github.com/AGSCourse)
[African Genomics Short Course Home Page at Wellcome Connecting Science](https://coursesandconferences.wellcomeconnectingscience.org/event/african-genomics-short-course-20251130/)
