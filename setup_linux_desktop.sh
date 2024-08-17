# This file sets up a Linux desktop entry, and moves plascad to the home directory.

printf "Moving the PlasCAD executable and icon to ~/plascad..."

chmod +x plascad

if [ ! -d ~/plascad ]; then
  mkdir ~/plascad
fi

cp plascad ~/plascad
cp icon.png ~/plascad/icon.png

# Update the desktop entry with the absolute path.
sed "s|~|$HOME|g" plascad.desktop > ~/.local/share/applications/plascad.desktop

printf "\nComplete. You can launch using plascad.desktop\n"