const fs = require("fs");
const path = require("path");

// Chromium aborts inside an AppImage when the OS blocks unprivileged user namespaces (a
// FUSE mount cannot hold a setuid helper), so a wrapper adds --no-sandbox only there.
exports.default = async function afterPack(context) {
    if (context.electronPlatformName !== "linux") return;

    const exe = context.packager.executableName;
    const exePath = path.join(context.appOutDir, exe);
    fs.renameSync(exePath, `${exePath}.bin`);

    const wrapper = `#!/bin/sh
# Only drop the Chromium sandbox where it cannot work at all: unprivileged user
# namespaces blocked (recent Ubuntu) and no setuid helper possible in an AppImage.
dir=$(dirname "$(readlink -f "$0")")
if [ "$(cat /proc/sys/kernel/apparmor_restrict_unprivileged_userns 2>/dev/null)" = "1" ] \\
   || [ "$(cat /proc/sys/kernel/unprivileged_userns_clone 2>/dev/null)" = "0" ]; then
    exec "$dir/${exe}.bin" --no-sandbox "$@"
fi
exec "$dir/${exe}.bin" "$@"
`;
    fs.writeFileSync(exePath, wrapper, {mode: 0o755});
};
