#!/bin/bash

# ======= 1. remove .git (remove all infos)=======
rm -rf .git

# ======= 2. Git User infos =======
git config --global user.name "li-wen-hua"
git config --global user.email "2990973334@qq.com"

# ======= 3. initialization =======
git init
git branch -m main # If no main branch exist

# ======= 4. add all files =======
git add .

# ======= 5. submmit to .git(Local Git) =======
git commit -m "Commit: 1. #Using SSH on $(date +%Y-%m-%d)"

# ======= 6. add remote Git address : SSH ("tokon" is forbbiden by School Internet) =======
git remote add origin git@github.com:li-wen-hua/LatticeReconstruct.git

# ======= 7. Push. Pull first and then Push, sometimes changes are made online =======
# git pull -u origin main --rebase
git push -u origin main

