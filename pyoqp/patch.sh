# apply a small patch to dl-find callback.py for dlf_get_multistate_gradients_wrapper

libdlfind_path=`python -m pip show libdlfind|grep Location|awk '{print $NF}'`
callback_file=$libdlfind_path/libdlfind/callback.py
patch_file='@@ -88,9 +88,11 @@\n         needcoupling: int,\n         iimage: int,\n         status: pointer[c_int],\n+        *args,\n+        **kwargs,\n     ) -> None:\n         coordinates_ = as_array(coords, (nvar,)).reshape(-1, 3)\n-        e_1, e_2, g_1, g_2, _ = func(coordinates_, needcoupling, iimage)\n+        e_1, e_2, g_1, g_2, _ = func(coordinates_, needcoupling, iimage, *args, **kwargs)\n         energy_ = as_array(energy, (2,))\n         energy_[0] = e_1\n         energy_[1] = e_2\n'

echo Found libdlfind $callback_file
echo Apply patch

printf "$patch_file" > libdlfind.patch
patch $callback_file libdlfind.patch
rm -rf libdlfind.patch

echo Done
